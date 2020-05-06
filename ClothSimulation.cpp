/***************************************************************************
# Copyright (c) 2015, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
***************************************************************************/
#include "glm/gtc/random.hpp"
#include "ClothSimulation.h"
#include "TriangleIntersection.h"

const char *ProjectTemplate::ClothVS = "Cloth.vs.hlsl";
const char *ProjectTemplate::ClothPS = "Cloth.ps.hlsl";
const char *ProjectTemplate::ModelPS = "Model.ps.hlsl";

const char *ProjectTemplate::SkyBoxTextures[] = {
    "Cubemaps/PiazzaDelPopolo.dds",
    "Cubemaps/SaintLazarusChurch.dds",
    "Cubemaps/SaintPetersSquare.dds",
    "Cubemaps/Sorsele3/Sorsele3.dds"

};

const char *ProjectTemplate::ModelUnitSphere = "Models/UnitSphere.fbx";
const char *ProjectTemplate::ModelUnitCylinder = "Models/UnitCylinder.fbx";
const char *ProjectTemplate::ModelUnitCone = "Models/UnitCone.fbx";

const int32_t ProjectTemplate::ClothSizeX = 32;
const int32_t ProjectTemplate::ClothSizeY = 32;

using float1 = float32_t;
using float3 = f32vec3;

void ClothPatch::init(int32_t width, int32_t height)
{
    mWidth = width;
    mHeight = height;
    mpParticles = new ClothParticle[width * height];

    for (int32_t y = 0; y < mHeight; ++y)
    {
        for (int32_t x = 0; x < mWidth; ++x)
        {
            auto *particle = getParticle(x, y);
            particle->mAcceleration = float3(0.0f, 0.0f, 0.0f);
            particle->mbStationary = false;
            particle->mPosition[0] = float32(x);
            particle->mPosition[1] = float32(y);
            particle->mPosition[2] = 0.0f;
            particle->mPrevPosition = particle->mPosition;
        }
    }
}

ClothParticle* ClothPatch::getParticle(int32_t x, int32_t y)
{
    return &mpParticles[y * mWidth + x];
}

CameraController& ProjectTemplate::getActiveCameraController()
{
    switch (mCameraType)
    {
    case ModelViewCamera:
        return mModelViewCameraController;
    case FirstPersonCamera:
        return mFirstPersonCameraController;
    case SixDoFCamera:
        return mSixDoFCameraController;
    default:
        should_not_get_here();
        return mModelViewCameraController;
    }
}

void ProjectTemplate::onGuiRender(SampleCallbacks* pSample, Gui* pGui)
{
    pGui->addText("Hello from ProjectTemplate");

    if (pGui->beginGroup("Visuals"))
    {

        Gui::DropdownList cameraDropdown;
        cameraDropdown.push_back({ ModelViewCamera, "Model-View" });
        cameraDropdown.push_back({ FirstPersonCamera, "First-Person" });
        cameraDropdown.push_back({ SixDoFCamera, "Six DoF" });
        pGui->addDropdown("Camera Type", cameraDropdown, reinterpret_cast<uint32_t&>(mCameraType));

        pGui->addCheckBox("Show Particles", mbShowParticles);
        pGui->addCheckBox("Show Wireframe", mbShowWireframe);
        pGui->addDirectionWidget("Light Direction", mLightDirection);
        pGui->endGroup();
    }
    pGui->addFloatVar("Air Temperature", mAirTemperature, -25.0f, 35.0f);
}

void ProjectTemplate::onLoad(SampleCallbacks* pSample, RenderContext* pRenderContext)
{
    srand(uint32_t(time(0)));

    const auto clothDiagonal = std::sqrtf(ClothSizeX * ClothSizeX + ClothSizeY * ClothSizeY);

    const auto vertexCount = ClothSizeX * ClothSizeY * 3;
    mpVBPositions = Buffer::create(vertexCount * sizeof(float32_t), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
    mpVBNormals = Buffer::create(vertexCount * sizeof(float32_t), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
    mpVBColors = Buffer::create(vertexCount * sizeof(float32_t), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);

    std::vector<uint16_t> vIndices;
    for (int32_t y = 0; y < ClothSizeY - 1; ++y)
    {
        for (int32_t x = 0; x < ClothSizeX - 1; ++x)
        {
            vIndices.push_back((y + 1) * ClothSizeY + x);
            vIndices.push_back(y * ClothSizeY + (x + 1));
            vIndices.push_back(y * ClothSizeY + x);

            vIndices.push_back((y + 1) * ClothSizeY + (x + 1));
            vIndices.push_back(y * ClothSizeY + (x + 1));
            vIndices.push_back((y + 1) * ClothSizeY + x);
        }
    }
    mpIndexBuffer = Buffer::create(vIndices.size() * sizeof(uint16_t), ResourceBindFlags::Index, Buffer::CpuAccess::None, vIndices.data());

    VertexLayout::SharedPtr pLayout = VertexLayout::create();
    VertexBufferLayout::SharedPtr pPositionsLayout = VertexBufferLayout::create();
    pPositionsLayout->addElement("POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
    pLayout->addBufferLayout(0, pPositionsLayout);
    VertexBufferLayout::SharedPtr pNormalsLayout = VertexBufferLayout::create();
    pNormalsLayout->addElement("NORMAL", 0, ResourceFormat::RGB32Float, 1, 1);
    pLayout->addBufferLayout(1, pNormalsLayout);
    VertexBufferLayout::SharedPtr pColorsLayout = VertexBufferLayout::create();
    pColorsLayout->addElement("COLOR", 0, ResourceFormat::RGB32Float, 1, 2);
    pLayout->addBufferLayout(2, pColorsLayout);

    Vao::BufferVec buffers{ mpVBPositions, mpVBNormals, mpVBColors };
    mpVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers, mpIndexBuffer, ResourceFormat::R16Uint);

    Program::DefineList defineList;

    GraphicsProgram::Desc clothProgramDesc;
    clothProgramDesc.addShaderLibrary(ClothVS).vsEntry("main");
    clothProgramDesc.addShaderLibrary(ClothPS).psEntry("main");
    mpClothProgram = GraphicsProgram::create(clothProgramDesc, defineList);
    mpClothVars = GraphicsVars::create(mpClothProgram->getReflector());
    mpClothState = GraphicsState::create();

    RasterizerState::Desc rasterizerDesc;
    rasterizerDesc.setCullMode(RasterizerState::CullMode::None);
    rasterizerDesc.setFillMode(RasterizerState::FillMode::Solid);
    mpRasterizeNormal = RasterizerState::create(rasterizerDesc);
    rasterizerDesc.setFillMode(RasterizerState::FillMode::Wireframe);
    mpRasterizeWireframe = RasterizerState::create(rasterizerDesc);

    mpClothState->setProgram(mpClothProgram);
    mpClothState->setVao(mpVao);

    mpModelProgram = GraphicsProgram::createFromFile(ModelPS, "", "main");
    mpModelVars = GraphicsVars::create(mpModelProgram->getReflector());
    mpModelState = GraphicsState::create();

    mpModelState->setProgram(mpModelProgram);

    mpCamera = Camera::create();

    mpCamera->setPosition(vec3(10, 1, 0));
    mpCamera->setTarget(vec3(0, 0, 0));
    mpCamera->setUpVector(vec3(0, 1, 0));
    mpCamera->setDepthRange(1.0f, 16.0f * 1024.0f);

    mpScene = Scene::create();
    mpSceneRenderer = SceneRenderer::create(mpScene);

    mLightDirection = vec3(0.13f, 0.27f, -0.9f);
    mpDirLight = DirectionalLight::create();
    mpDirLight->setWorldDirection(mLightDirection);

    mModelViewCameraController.attachCamera(mpCamera);
    mFirstPersonCameraController.attachCamera(mpCamera);
    mSixDoFCameraController.attachCamera(mpCamera);

    mModelViewCameraController.setModelParams(vec3(0.0f), clothDiagonal / 2.0f, 2.0f);
    mFirstPersonCameraController.setCameraSpeed(clothDiagonal * 0.25f);
    mSixDoFCameraController.setCameraSpeed(clothDiagonal * 0.25f);

    mCameraType = ModelViewCamera;

    mPrevTime = pSample->getCurrentTime();
    mClothPatch.init(ClothSizeX, ClothSizeY);
    mHeldParticle = nullptr;

    mClothPatch.getParticle(0, 0)->mbStationary = true;
    mClothPatch.getParticle(ClothSizeX - 1, 0)->mbStationary = true;

    mpUnitSphere = Model::createFromFile(ModelUnitSphere, Model::LoadFlags::None);
    mpUnitCylinder = Model::createFromFile(ModelUnitCylinder, Model::LoadFlags::None);
    mpUnitCone = Model::createFromFile(ModelUnitCone, Model::LoadFlags::None);

    mpParticleDbgModels = new ParticleDbgModels[mClothPatch.mWidth * mClothPatch.mHeight];
    for (int32_t y = 0; y < mClothPatch.mHeight; ++y)
    {
        for (int32_t x = 0; x < mClothPatch.mWidth; ++x)
        {
            auto pParticle = mClothPatch.getParticle(x, y);
            auto &dbgModels = mpParticleDbgModels[y * mClothPatch.mWidth + x];

            dbgModels.pUnitSphere = Scene::ModelInstance::create(mpUnitSphere, pParticle->mPosition, vec3(0.0f), vec3(0.001f));
            dbgModels.pUnitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
            mpScene->addModelInstance(dbgModels.pUnitSphere);

            dbgModels.pUnitCylinder = Scene::ModelInstance::create(mpUnitCylinder, pParticle->mPosition, vec3(0.0f), vec3(0.1f));
            dbgModels.pUnitCylinder->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
            mpScene->addModelInstance(dbgModels.pUnitCylinder);
        }
    }

    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(Sampler::Filter::Linear, Sampler::Filter::Linear, Sampler::Filter::Linear);
    mpTriLinearSampler = Sampler::create(samplerDesc);

    mpSkybox = SkyBox::create(SkyBoxTextures[linearRand<int32>(0, arraysize(SkyBoxTextures) - 1)], true, mpTriLinearSampler);
}

static constexpr struct
{
    float32 mTemperature, mDensity;
} AirPropertiesTable[] = {
    { -25.0f, 1.4224f },
    { -20.0f, 1.3943f },
    { -15.0f, 1.3673f },
    { -10.0f, 1.3413f },
    { - 5.0f, 1.3163f },
    {   0.0f, 1.2922f },
    {   5.0f, 1.2690f },
    {  10.0f, 1.2466f },
    {  15.0f, 1.2250f },
    {  20.0f, 1.2041f },
    {  25.0f, 1.1839f },
    {  30.0f, 1.1644f },
    {  35.0f, 1.1455f },
};

static float32 GetAirDensity(float32 temperature)
{
    temperature = glm::clamp(
        temperature,
        AirPropertiesTable[0].mTemperature,
        AirPropertiesTable[arraysize(AirPropertiesTable) - 1].mTemperature
    );

    for (int32_t i = 1; i < arraysize(AirPropertiesTable); ++i)
    {
        const auto &temperature0 = AirPropertiesTable[i - 1].mTemperature;
        const auto &temperature1 = AirPropertiesTable[i].mTemperature;

        if (temperature <= temperature1)
        {
            return lerp(
                AirPropertiesTable[i - 1].mDensity,
                AirPropertiesTable[i].mDensity,
                (temperature - temperature0) / (temperature1 - temperature0)
            );
        }
    }

    return AirPropertiesTable[0].mDensity;
}

static vec3 getYawPitchRoll(vec3 direction)
{
    direction = normalize(direction);
    const auto phi = atan2(direction[0], direction[2]);
    const auto theta = asin(-direction[1]);
    return vec3(phi, theta + radians(90.0f), 0.0f);
}

void ProjectTemplate::onFrameRender(SampleCallbacks* pSample, RenderContext* pRenderContext, const Fbo::SharedPtr& pTargetFbo)
{
    const vec4 clearColor(0.38f, 0.52f, 0.10f, 1);
    pRenderContext->clearFbo(pTargetFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::All);

    mpSkybox->render(pRenderContext, mpCamera.get());

    const auto dragCoefficient = 0.47f;
    const auto area = 1.0f;
    const auto airDensity = GetAirDensity(mAirTemperature);

    const auto windAcceleration = float3(3.0f, 0.0f, 3.0f) / mClothPatch.mMass;
    const auto gravityAcceleration = float3(0.0f, -9.8067f, 0.0f);
    const auto currTime = pSample->getCurrentTime();
    const auto deltaTime = currTime - mPrevTime;
    mPrevTime = currTime;

    for (int32_t y = 0; y < mClothPatch.mHeight; ++y)
    {
        for (int32_t x = 0; x < mClothPatch.mWidth; ++x)
        {
            auto particle = mClothPatch.getParticle(x, y);
            if (particle->mbStationary) continue;

            const auto currPosition = particle->mPosition;
            const auto velocity = particle->mPosition - particle->mPrevPosition;

            const auto dragForce = 0.5f * airDensity * dragCoefficient * area * (velocity * velocity);
            const auto dragAcceleration = dragForce * mClothPatch.mMass;

            const auto acceleration = gravityAcceleration + windAcceleration /*+ dragAcceleration*/;

            particle->mPosition += velocity + acceleration * deltaTime * deltaTime;
            particle->mPrevPosition = currPosition;

            if (mbShowParticles)
            {
                auto &dbgModels = mpParticleDbgModels[y * mClothPatch.mWidth + x];
                dbgModels.pUnitSphere->setTranslation(particle->mPosition, true);
                dbgModels.pUnitCylinder->setTranslation(particle->mPosition, true);
                dbgModels.pUnitCylinder->setRotation(getYawPitchRoll(velocity));
                dbgModels.pUnitCylinder->setScaling(vec3(0.1, 0.1, length(velocity)));
            }
        }
    }

    const int32_t Offsets[][2] = {
        {  0,  1 },
        {  0, -1 },
        {  1,  0 },
        { -1,  0 },

        {  1,  1 },
        {  1, -1 },
        { -1,  1 },
        { -1, -1 },

        {  0,  2 },
        {  0, -2 },
        {  2,  0 },
        { -2,  0 },
    };

    for (int32_t y = 0; y < mClothPatch.mHeight; ++y)
    {
        for (int32_t x = 0; x < mClothPatch.mWidth; ++x)
        {
            auto first = mClothPatch.getParticle(x, y);
            float3 adjust(0.0f);

            for (const auto &offset : Offsets)
            {
                if (x + offset[0] < 0 || x + offset[0] >= mClothPatch.mWidth) continue;
                if (y + offset[1] < 0 || y + offset[1] >= mClothPatch.mHeight) continue;

                auto second = mClothPatch.getParticle(x + offset[0], y + offset[1]);
                auto distance = first->mPosition - second->mPosition;
                auto current = length(distance);
                auto target = sqrtf(float32_t(offset[0] * offset[0] + offset[1] * offset[1]));

                float weightFirst, weightSecond;
                if (first->mbStationary && second->mbStationary)
                    weightFirst = 0.0f, weightSecond = 0.0f;
                else if (!first->mbStationary && second->mbStationary)
                    weightFirst = 1.0f, weightSecond = 0.0f;
                else if (first->mbStationary && !second->mbStationary)
                    weightFirst = 0.0f, weightSecond = 1.0f;
                else
                    weightFirst = 0.5f, weightSecond = 0.5f;

                distance *= ((target - current) / current);
                first->mPosition += distance * weightFirst;
                second->mPosition -= distance * weightSecond;
            }
        }
    }

    float32_t *pPositions = reinterpret_cast<float32_t*>(mpVBPositions->map(Buffer::MapType::WriteDiscard));
    float32_t *pNormals = reinterpret_cast<float32_t*>(mpVBNormals->map(Buffer::MapType::WriteDiscard));
    for (int32_t y = 0; y < mClothPatch.mHeight; ++y)
    {
        for (int32_t x = 0; x < mClothPatch.mWidth; ++x)
        {
            const auto *particle = mClothPatch.getParticle(x, y);
            auto *position = pPositions + (y * mClothPatch.mWidth + x) * 3;

            position[0] = particle->mPosition[0];
            position[1] = particle->mPosition[1];
            position[2] = particle->mPosition[2];

            // TODO: Fit a plane to the particles' position and use its normal
            
            // Normal sampling pattern:
            // *  *  *
            //   /|\
            //  /4|3\
            // *__*__*
            //  \1|2/
            //   \|/
            // *  *  *

            const int32_t NormalSamples[][3][2] =
            {
                { { -1,  0 }, {  0,  0 }, {  0,  1 } }, // 1
                { {  0,  1 }, {  0,  0 }, {  1,  0 } }, // 2
                { {  1,  0 }, {  0,  0 }, {  0, -1 } }, // 3
                { {  0, -1 }, {  0,  0 }, { -1,  0 } }, // 4
            };

            const auto &GetParticleSample = [&](const int32_t(&offset)[2]) -> float3 {

                int32_t x1 = glm::clamp(x + offset[0], 0, mClothPatch.mWidth);
                int32_t y1 = glm::clamp(y + offset[1], 0, mClothPatch.mHeight);
                return mClothPatch.getParticle(x1, y1)->mPosition;
            };
            const auto &GetPlaneNormal = [](const float3 &A, const float3 &B, const float3 &C) -> float3 {
                const auto result = normalize(cross(C - B, A - B));
                return any(isnan(result)) ? float3(0.0f) : result;
            };

            float3 accumulator(0);
            for (const auto &sample : NormalSamples)
            {
                accumulator += GetPlaneNormal(
                    GetParticleSample(sample[0]), GetParticleSample(sample[1]), GetParticleSample(sample[2])
                );
            }

            accumulator = normalize(accumulator);
            auto outNormal = pNormals + (y * mClothPatch.mWidth + x) * 3;
            outNormal[0] = accumulator[0];
            outNormal[1] = accumulator[1];
            outNormal[2] = accumulator[2];
        }
    }
    mpVBPositions->unmap();
    mpVBNormals->unmap();

    const auto &GetTriangles = [&](int32_t x, int32_t y, vec3 triangle0[3], vec3 triangle1[3]) -> void
    {
        const auto &GetVertex = [&](int32_t x, int32_t y) -> vec3
        {
            return mClothPatch.getParticle(x, y)->mPosition;
        };

        triangle0[0] = GetVertex(x, y + 1);
        triangle0[1] = GetVertex(x + 1, y);
        triangle0[2] = GetVertex(x, y);

        triangle1[0] = GetVertex(x + 1, y + 1);
        triangle1[1] = GetVertex(x + 1, y);
        triangle1[2] = GetVertex(x, y + 1);
    };

    auto *pColors = reinterpret_cast<vec3*>(mpVBColors->map(Buffer::MapType::WriteDiscard));
    std::fill_n(pColors, mClothPatch.mHeight * mClothPatch.mWidth, vec3(0.0f, 1.0f, 0.0f));

    const auto &HandleCollision = [&](int32_t x, int32_t y) -> void
    {
        pColors[y * mClothPatch.mWidth + x] = vec3(1.0f, 0.0f, 0.0f);
        // mClothPatch.getParticle(x, y)->mPosition = mClothPatch.getParticle(x, y)->mPrevPosition;
    };

    for (int32_t y0 = 0; y0 < mClothPatch.mHeight - 1; ++y0)
    {
        for (int32_t x0 = 0; x0 < mClothPatch.mWidth - 1; ++x0)
        {
            vec3 triangle0_0[3], triangle0_1[3];
            GetTriangles(x0, y0, triangle0_0, triangle0_1);

            for (int32_t y1 = y0; y1 < mClothPatch.mHeight - 1; ++y1)
            {
                for (int32_t x1 = y1 == y0 ? x0 + 1 : 0; x1 < mClothPatch.mWidth - 1; ++x1)
                {
                    if (std::abs(y0 - y1) <= 1 && std::abs(x0 - x1) <= 1)
                        continue;

                    vec3 triangle1_0[3], triangle1_1[3];
                    GetTriangles(x1, y1, triangle1_0, triangle1_1);

                    vec3 responses[4];

                    bool collisions[4] =
                    {
                        intersectTriangles(triangle0_0, triangle1_0, responses[0]),
                        intersectTriangles(triangle0_0, triangle1_1, responses[1]),
                        intersectTriangles(triangle0_1, triangle1_0, responses[2]),
                        intersectTriangles(triangle0_1, triangle1_1, responses[3])
                    };

                    if (collisions[0] || collisions[1])
                    {
                        HandleCollision(x0, y0 + 1);
                        HandleCollision(x0 + 1, y0);
                        HandleCollision(x0, y0);
                    }
                    if (collisions[2] || collisions[3])
                    {
                        HandleCollision(x0 + 1, y0 + 1);
                        HandleCollision(x0 + 1, y0);
                        HandleCollision(x0, y0 + 1);
                    }
                    if (collisions[0] || collisions[2])
                    {
                        HandleCollision(x1, y1 + 1);
                        HandleCollision(x1 + 1, y1);
                        HandleCollision(x1, y1);
                    }
                    if (collisions[1] || collisions[3])
                    {
                        HandleCollision(x1 + 1, y1 + 1);
                        HandleCollision(x1 + 1, y1);
                        HandleCollision(x1, y1 + 1);
                    }
                }
            }
        }
    }
    mpVBColors->unmap();

    const auto showVectorField = [&](BoundingBox box, vec3 step, vec3 direction) -> void
    {
        auto scene = Scene::create();
        auto sceneRenderer = SceneRenderer::create(scene);

        const auto min = box.getMinPos();
        const auto max = box.getMaxPos();
        const auto yawPitchRoll = getYawPitchRoll(direction);

        for (auto x = min.x; x <= max.x; x += step.x)
        {
            for (auto y = min.y; y <= max.y; y += step.y)
            {
                for (auto z = min.z; z <= max.z; z += step.z)
                {
                    auto unitCylinder = Scene::ModelInstance::create(mpUnitCylinder, vec3(x, y, z), yawPitchRoll, vec3(0.1f));
                    unitCylinder->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
                    scene->addModelInstance(unitCylinder);

                    auto unitCone = Scene::ModelInstance::create(mpUnitCone, vec3(x, y, z) + direction, yawPitchRoll, vec3(0.1f));
                    unitCone->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
                    scene->addModelInstance(unitCone);
                }
            }
        }

        sceneRenderer->update(currTime);
        mpDirLight->setIntoProgramVars(mpModelVars.get(), mpModelVars["PerFrameCB"].get(), "gDirLight");

        pRenderContext->setGraphicsVars(mpModelVars);
        pRenderContext->setGraphicsState(mpModelState);

        mpModelState->pushFbo(pTargetFbo);
        sceneRenderer->renderScene(pRenderContext, mpCamera.get());
        mpModelState->popFbo();
    };

    if (mbShowWireframe)
        mpClothState->setRasterizerState(mpRasterizeWireframe);
    else
        mpClothState->setRasterizerState(mpRasterizeNormal);
    mpDirLight->setWorldDirection(mLightDirection);

    getActiveCameraController().update();
    mpCamera->setIntoConstantBuffer(mpClothVars->getConstantBuffer("InternalPerFrameCB").get(), "gCamera");
    mpDirLight->setIntoProgramVars(mpClothVars.get(), mpClothVars["PerFrameCB"].get(), "gDirLight");

    pRenderContext->setGraphicsVars(mpClothVars);
    pRenderContext->setGraphicsState(mpClothState);

    mpClothState->pushFbo(pTargetFbo);
    pRenderContext->drawIndexed((mClothPatch.mWidth - 1) * (mClothPatch.mHeight - 1) * 6, 0, 0);
    mpClothState->popFbo();

    //BoundingBox box;
    //box.center = vec3(0, 0, 0);
    //box.extent = vec3(40, 40, 40);
    //showVectorField(box, vec3(10, 10, 10), mLightDirection);

    if (mbShowParticles)
    {
        mpSceneRenderer->update(currTime);
        mpDirLight->setIntoProgramVars(mpModelVars.get(), mpModelVars["PerFrameCB"].get(), "gDirLight");

        pRenderContext->setGraphicsVars(mpModelVars);
        pRenderContext->setGraphicsState(mpModelState);

        mpModelState->pushFbo(pTargetFbo);
        mpSceneRenderer->renderScene(pRenderContext, mpCamera.get());
        mpModelState->popFbo();
    }
}

void ProjectTemplate::onShutdown(SampleCallbacks* pSample)
{
}

bool ProjectTemplate::onKeyEvent(SampleCallbacks* pSample, const KeyboardEvent& keyEvent)
{
    return getActiveCameraController().onKeyEvent(keyEvent);
}

bool ProjectTemplate::onMouseEvent(SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
    const auto &cameraData = mpCamera->getData();
    const auto &GetMouseRay = [&](vec3 &rayOrigin, vec3 &rayDirection) -> void
    {
        rayOrigin = cameraData.posW;
        rayDirection = cameraData.cameraW;
        rayDirection += cameraData.cameraU * (mouseEvent.pos.x - 0.5f) * 2.0f;
        rayDirection += cameraData.cameraV * (0.5f - mouseEvent.pos.y) * 2.0f;
        rayDirection = normalize(rayDirection);
    };

    switch (mouseEvent.type)
    {
    case MouseEvent::Type::LeftButtonDown:
    {
        if (mHeldParticle == nullptr)
        {
            vec3 rayOrigin, rayDirection;
            GetMouseRay(rayOrigin, rayDirection);

            const auto &GetTriangles = [&](int32_t x, int32_t y, vec3 triangle0[3], vec3 triangle1[3]) -> void
            {
                const auto &GetVertex = [&](int32_t x, int32_t y) -> vec3
                {
                    return mClothPatch.getParticle(x, y)->mPosition;
                };

                triangle0[0] = GetVertex(x, y + 1);
                triangle0[1] = GetVertex(x + 1, y);
                triangle0[2] = GetVertex(x, y);

                triangle1[0] = GetVertex(x + 1, y + 1);
                triangle1[1] = GetVertex(x + 1, y);
                triangle1[2] = GetVertex(x, y + 1);
            };

            auto closestHit = std::numeric_limits<float32>::max();

            for (int32_t y = 0; y < mClothPatch.mHeight - 1; ++y)
            {
                for (int32_t x = 0; x < mClothPatch.mWidth - 1; ++x)
                {
                    vec3 triangle0[3], triangle1[3], intersection;
                    GetTriangles(x, y, triangle0, triangle1);
                    auto currentHit = 0.0f;

                    if (intersectTriangle(triangle0, rayOrigin, rayDirection, intersection))
                    {
                        if ((currentHit = length(intersection - triangle0[0])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x, y + 1);
                            currentHit = closestHit;
                        }
                        if ((currentHit = length(intersection - triangle0[1])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x + 1, y);
                            currentHit = closestHit;
                        }
                        if ((currentHit = length(intersection - triangle0[2])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x, y);
                            currentHit = closestHit;
                        }
                    }
                    if (intersectTriangle(triangle1, rayOrigin, rayDirection, intersection))
                    {
                        if ((currentHit = length(intersection - triangle1[0])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x + 1, y + 1);
                            currentHit = closestHit;
                        }
                        if ((currentHit = length(intersection - triangle1[1])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x + 1, y);
                            currentHit = closestHit;
                        }
                        if ((currentHit = length(intersection - triangle1[2])) < closestHit)
                        {
                            mHeldParticle = mClothPatch.getParticle(x, y + 1);
                            currentHit = closestHit;
                        }
                    }
                }
            }

            if (mHeldParticle != nullptr)
            {
                mHeldParticle->mbStationary = true;
                return false;
            }
        }
    } break;
    case MouseEvent::Type::LeftButtonUp:
    {
        if (mHeldParticle != nullptr)
        {
            mHeldParticle->mbStationary = false;
            mHeldParticle = nullptr;
            return false;
        }

    } break;
    }

    if (mHeldParticle != nullptr)
    {
        vec3 rayOrigin, rayDirection;
        GetMouseRay(rayOrigin, rayDirection);

        const auto UxV = cross(cameraData.cameraU, cameraData.cameraV);
        const auto t = dot(rayOrigin - mHeldParticle->mPosition, UxV) / dot(-rayDirection, UxV);
        mHeldParticle->mPosition = rayOrigin + t * rayDirection;
    }

    return getActiveCameraController().onMouseEvent(mouseEvent);
}

void ProjectTemplate::onDataReload(SampleCallbacks* pSample)
{

}

void ProjectTemplate::onResizeSwapChain(SampleCallbacks* pSample, uint32_t width, uint32_t height)
{
    mpCamera->setFocalLength(21.0f);
    float32_t aspectRatio = float32_t(width) / float32_t(height);
    mpCamera->setAspectRatio(aspectRatio);
}

int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nShowCmd)
{
    ProjectTemplate::UniquePtr pRenderer = std::make_unique<ProjectTemplate>();
    SampleConfig config;
    config.windowDesc.title = "Falcor Project Template";
    config.windowDesc.resizableWindow = true;
    Sample::run(config, pRenderer);
    return 0;
}
