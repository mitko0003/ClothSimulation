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
#include "ClothSample.h"

const char *ClothSample::ModelPS = "Model.ps.hlsl";

const char *ClothSample::SkyBoxTextures[] = {
    "Cubemaps/PiazzaDelPopolo.dds",
    "Cubemaps/SaintLazarusChurch.dds",
    "Cubemaps/SaintPetersSquare.dds",
    "Cubemaps/Sorsele3/Sorsele3.dds"
};

const int32_t ClothSample::ClothSizeX = 32;
const int32_t ClothSample::ClothSizeY = 32;

CameraController& ClothSample::getActiveCameraController()
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

void ClothSample::onGuiRender(SampleCallbacks* pSample, Gui* pGui)
{
    pGui->addText("Hello from ProjectTemplate");

    if (pGui->beginGroup("Visuals"))
    {

        Gui::DropdownList cameraDropdown;
        cameraDropdown.push_back({ ModelViewCamera, "Model-View" });
        cameraDropdown.push_back({ FirstPersonCamera, "First-Person" });
        cameraDropdown.push_back({ SixDoFCamera, "Six DoF" });
        pGui->addDropdown("Camera Type", cameraDropdown, reinterpret_cast<uint32_t&>(mCameraType));

        pGui->addDirectionWidget("Light Direction", mLightDirection);
        pGui->endGroup();
    }
	mClothPatch.onGuiRender(this, pSample);
    pGui->addFloatVar("Air Temperature", mAirTemperature, -25.0f, 35.0f);
}

void ClothSample::onLoad(SampleCallbacks* pSample, RenderContext* pRenderContext)
{
    srand(uint32_t(time(0)));

    mpModelProgram = GraphicsProgram::createFromFile(ModelPS, "", "main");
    mpModelVars = GraphicsVars::create(mpModelProgram->getReflector());
    mpModelState = GraphicsState::create();

    mpModelState->setProgram(mpModelProgram);

    mpCamera = Camera::create();

    mpCamera->setPosition(vec3(10, 1, 0));
    mpCamera->setTarget(vec3(0, 0, 0));
    mpCamera->setUpVector(vec3(0, 1, 0));
    mpCamera->setDepthRange(1.0f, 16.0f * 1024.0f);

    mLightDirection = vec3(0.13f, 0.27f, -0.9f);
    mpDirLight = DirectionalLight::create();
    mpDirLight->setWorldDirection(mLightDirection);

    mModelViewCameraController.attachCamera(mpCamera);
    mFirstPersonCameraController.attachCamera(mpCamera);
    mSixDoFCameraController.attachCamera(mpCamera);

	const auto clothDiagonal = 50.0f; // update this based on mesh?
    mModelViewCameraController.setModelParams(vec3(0.0f), clothDiagonal / 2.0f, 2.0f);
    mFirstPersonCameraController.setCameraSpeed(clothDiagonal * 0.25f);
    mSixDoFCameraController.setCameraSpeed(clothDiagonal * 0.25f);

    mCameraType = ModelViewCamera;

    mPrevTime = pSample->getCurrentTime();
    mClothPatch.init(ClothSizeX, ClothSizeY);

    mClothPatch.getParticle(0, 0)->mbStationary = true;
    mClothPatch.getParticle(ClothSizeX - 1, 0)->mbStationary = true;

    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(Sampler::Filter::Linear, Sampler::Filter::Linear, Sampler::Filter::Linear);
    mpTriLinearSampler = Sampler::create(samplerDesc);

    mpSkybox = SkyBox::create(SkyBoxTextures[linearRand<int32>(0, arraysize(SkyBoxTextures) - 1)], true, mpTriLinearSampler);
}

void ClothSample::onFrameRender(SampleCallbacks* pSample, RenderContext* pRenderContext, const Fbo::SharedPtr& pTargetFbo)
{
    const vec4 clearColor(0.38f, 0.52f, 0.10f, 1);
    pRenderContext->clearFbo(pTargetFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::All);

	const auto currTime = pSample->getCurrentTime();
	const auto deltaTime = currTime - mPrevTime;
	mPrevTime = currTime;

	mClothPatch.simulate(deltaTime);
    mpSkybox->render(pRenderContext, mpCamera.get());

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

    mpDirLight->setWorldDirection(mLightDirection);

    getActiveCameraController().update();

    //BoundingBox box;
    //box.center = vec3(0, 0, 0);
    //box.extent = vec3(40, 40, 40);
    //showVectorField(box, vec3(10, 10, 10), mLightDirection);
}

void ClothSample::onShutdown(SampleCallbacks* pSample)
{
}

bool ClothSample::onKeyEvent(SampleCallbacks* pSample, const KeyboardEvent& keyEvent)
{
    return getActiveCameraController().onKeyEvent(keyEvent);
}

bool ClothSample::onMouseEvent(SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
	if (!mClothPatch.onMouseEvent(this, pSample, mouseEvent))
		return getActiveCameraController().onMouseEvent(mouseEvent);
	return false;
}

void ClothSample::onDataReload(SampleCallbacks* pSample)
{

}

void ClothSample::onResizeSwapChain(SampleCallbacks* pSample, uint32_t width, uint32_t height)
{
    mpCamera->setFocalLength(21.0f);
    float32_t aspectRatio = float32_t(width) / float32_t(height);
    mpCamera->setAspectRatio(aspectRatio);
}

int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nShowCmd)
{
    ClothSample::UniquePtr pRenderer = std::make_unique<ClothSample>();
    SampleConfig config;
    config.windowDesc.title = "Falcor Project Template";
    config.windowDesc.resizableWindow = true;
    Sample::run(config, pRenderer);
    return 0;
}
