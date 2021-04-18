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
#include "ClothSample.h"
#include "TriangleIntersection.h"
#include "glm/gtc/random.hpp"
#include "glm/gtx/projection.hpp"
#include "Externals/dear_imgui/imgui.h"

const char *ClothSample::ModelPS = "Model.ps.hlsl";

const char *ClothSample::ModelUnitSphere = "Models/UnitSphere.fbx";
const char *ClothSample::ModelUnitCylinder = "Models/UnitCylinder.fbx";
const char *ClothSample::ModelUnitCone = "Models/UnitCone.fbx";
#pragma optimize("", off)
const char *ClothSample::SkyBoxTextures[] = {
    "Cubemaps/PiazzaDelPopolo.dds",
    "Cubemaps/SaintLazarusChurch.dds",
    "Cubemaps/SaintPetersSquare.dds",
    "Cubemaps/Sorsele3/Sorsele3.dds"
};

void PhysicsObject::testSelection(SelectionQuery &query)
{
    auto distance = intersect(query.ray);
    if (query.closestHit > distance)
    {
        query.closestHit = distance;
        query.closestObject = this;
    }
}

void PhysicsObject::makeSelection(ClothSample *pClothSample, SelectionQuery&)
{
    mpController = ObjectController::SharedPtr(new ObjectController(pClothSample, getPosition(), getScale()));
    pClothSample->mpObjectController = mpController;
}

void PhysicsObject::loseSelection(ClothSample *pClothSample)
{
    mpController = nullptr;
    pClothSample->mpObjectController = nullptr;
}

bool PhysicsObject::onMouseEvent(ClothSample *pClothSample, SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
    if (mpController->onMouseEvent(pClothSample, pSample, mouseEvent))
    {
        applyController(mpController.get());
        return true;
    }
    return false;
}

class PhysicsSphere final : public PhysicsObject
{
public:
    PhysicsSphere(glm::vec3 center, float radius) :
        center(center), radius(radius) {}

    vec3 getPosition() const override
    {
        return center;
    }

    vec3 getScale() const override
    {
        return vec3(radius);
    }

    // Courtesy of: https://fiftylinesofcode.com/ray-sphere-intersection/
    float intersect(const Ray &ray) override
    {
        vec3 originToCenter = ray.origin - center;

        float p = glm::dot(ray.direction, originToCenter);
        float q = glm::dot(originToCenter, originToCenter) - (radius * radius);

        float discriminant = (p * p) - q;
        if (discriminant < 0.0f)
            return std::numeric_limits<float>::infinity();

        float dRoot = sqrt(discriminant);
        return std::min(-p - dRoot, -p + dRoot);
    }

    vec3 intersectionResponse(vec3 V[3]) override
    {
        const auto closestPoint = triangleClosestPoint(V, center);
        const auto diff = distance(closestPoint, center);
        if (diff < radius)
        {
            const auto direction = normalize(closestPoint - center);
            return direction * (radius - diff);
        }
        return vec3(0.0f);
    }

    void render(ClothSample *pClothSample, SampleCallbacks*) override
    {
        auto unitSphere = Scene::ModelInstance::create(pClothSample->mpDbgUnitSphere, center, vec3(0.0f), vec3(radius * 0.01f));
        unitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
        pClothSample->mpScene->addModelInstance(unitSphere);
    }

    void applyController(const ObjectController *controller) override
    {
        center = controller->getPosition();
        radius = std::max(controller->getScale().x, 0.0001f);
    }

private:
    vec3 center;
    float radius;
};

glm::vec3 ObjectController::getPosition() const
{
    return mpProxyModelInstance->getTranslation();
}

glm::vec3 ObjectController::getScale() const
{
    return mpProxyModelInstance->getScaling();
}

ObjectController::ObjectController(ClothSample *pClothSample, const glm::vec3 &translation, const glm::vec3 &scale)
{
    auto backBufferFBO = gpDevice->getSwapChainFbo();
    uint32_t backBufferWidth = backBufferFBO->getWidth();
    uint32_t backBufferHeight = backBufferFBO->getHeight();

    mpScene = Scene::create();
    mpSceneRenderer = SceneRenderer::create(mpScene);
    mpPicker = Picking::create(mpScene, backBufferWidth, backBufferHeight);

    mpProxyModelInstance = Scene::ModelInstance::create(pClothSample->mpDbgUnitSphere, translation, vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f), scale);
    mpProxyModelInstance->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
    mpScene->addModelInstance(mpProxyModelInstance);

    mGizmos[int32_t(Gizmo::Type::Translate)] = ::TranslateGizmo::create(mpScene, "Framework/Models/TranslateGizmo.obj");
    mGizmos[int32_t(Gizmo::Type::Rotate)] = ::RotateGizmo::create(mpScene, "Framework/Models/RotateGizmo.obj");
    mGizmos[int32_t(Gizmo::Type::Scale)] = ::ScaleGizmo::create(mpScene, "Framework/Models/ScaleGizmo.obj");

    const auto &activeGizmo = mGizmos[int32_t(mActiveGizmoType)];
    activeGizmo->setVisible(true);
    activeGizmo->setTransform(pClothSample->mpCamera, mpProxyModelInstance);
}

void ObjectController::renderScene(ClothSample *pClothSample, SampleCallbacks *pSample)
{
    const auto currTime = pSample->getCurrentTime();
    auto *pRenderContext = pSample->getRenderContext();

    const auto &activeGizmo = mGizmos[int32_t(mActiveGizmoType)];
    activeGizmo->setTransform(pClothSample->mpCamera, mpProxyModelInstance);

    const auto cashedScale = mpProxyModelInstance->getScaling();
    mpProxyModelInstance->setScaling(vec3(0.0001f));

    mpSceneRenderer->update(currTime);
    mpSceneRenderer->renderScene(pRenderContext, pClothSample->mpCamera.get());

    mpProxyModelInstance->setScaling(cashedScale);
}

void ObjectController::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
    auto *pGui = pSample->getGui();

    if (pGui->beginGroup("Object Controller", true))
    {
        Gui::DropdownList gizmoTypeDropdown;
        gizmoTypeDropdown.push_back({ int32_t(Gizmo::Type::Translate), "Translate" });
        gizmoTypeDropdown.push_back({ int32_t(Gizmo::Type::Rotate), "Rotate" });
        gizmoTypeDropdown.push_back({ int32_t(Gizmo::Type::Scale), "Scale" });

        const auto &oldGizmo = mGizmos[int32_t(mActiveGizmoType)];
        if (pGui->addDropdown("Gizmo Type", gizmoTypeDropdown, reinterpret_cast<uint32_t&>(mActiveGizmoType)))
        {
            oldGizmo->setVisible(false);
            const auto &newGizmo = mGizmos[int32_t(mActiveGizmoType)];
            newGizmo->setVisible(true);
            newGizmo->setTransform(pClothSample->mpCamera, mpProxyModelInstance);
        }

        pGui->endGroup();
    }
}

bool ObjectController::onMouseEvent(ClothSample *pClothSample, SampleCallbacks *pSample, const MouseEvent &mouseEvent)
{
    switch (mouseEvent.type)
    {
    case MouseEvent::Type::LeftButtonDown:
        // Gizmo Selection
        if (mGizmoBeingDragged == false)
        {
            const auto &pCamera = pClothSample->mpCamera;
            auto *pRenderContext = pSample->getRenderContext();

            if (mpPicker->pick(pRenderContext, mouseEvent.pos, pCamera))
            {
                const auto& pInstance = mpPicker->getPickedModelInstance();

                // If picked model instance is part of the active gizmo
                if (mGizmos[int32_t(mActiveGizmoType)]->beginAction(pCamera, pInstance))
                {
                    mGizmoBeingDragged = true;
                    mGizmos[int32_t(mActiveGizmoType)]->update(pCamera, mouseEvent);
                    return true;
                }
            }
        }
        break;

    case MouseEvent::Type::Move:
        // Gizmo Drag
        if (mGizmoBeingDragged)
        {
            const auto &pCamera = pClothSample->mpCamera;
            mGizmos[(uint32_t)mActiveGizmoType]->update(pCamera, mouseEvent);

            const auto &activeGizmo = mGizmos[int32_t(mActiveGizmoType)];
            activeGizmo->applyDelta(mpProxyModelInstance);

            if (mActiveGizmoType == Gizmo::Type::Scale && mUniformScale)
            {
                auto scale = mpProxyModelInstance->getScaling();
                     if (scale.x == scale.y) scale.x = scale.y = scale.z;
                else if (scale.x == scale.z) scale.x = scale.z = scale.y;
                else if (scale.y == scale.z) scale.y = scale.z = scale.x;
                mpProxyModelInstance->setScaling(scale);
            }

            return true;
        }
        break;

    case MouseEvent::Type::LeftButtonUp:
        if (mGizmoBeingDragged)
        {
            mGizmoBeingDragged = false;
            return true;
        }
        break;
    }

    return false;
}

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

void ClothSample::updateClothModel()
{
	mpClothModel = ClothModel::createClothModel(mClothModelType, mpClothModel.get());
	mpClothModel->init(mClothPatchSize, mClothPatchTessellation);
}

void ClothSample::onGuiRender(SampleCallbacks* pSample, Gui* pGui)
{
    if (pGui->beginGroup("Visuals"))
    {
        Gui::DropdownList cameraDropdown;
        cameraDropdown.push_back({ ModelViewCamera, "Model-View" });
        cameraDropdown.push_back({ FirstPersonCamera, "First-Person" });
        cameraDropdown.push_back({ SixDoFCamera, "Six DoF" });
        pGui->addDropdown("Camera Type", cameraDropdown, reinterpret_cast<uint32_t&>(mCameraType));

		if (pGui->addFloat3Slider("Light Direction", mLightDirection, -1.0f, 1.0f))
		{
			if (length(mLightDirection) != 0.0f)
				mLightDirection = normalize(mLightDirection);
			showDirectionHelper(mLightDirection);
		}

        pGui->endGroup();
    }

    if (mpObjectController != nullptr)
        mpObjectController->onGuiRender(this, pSample);

	if (pGui->beginGroup("Cloth Setup", true))
	{
        if (pGui->addButton("Restart"))
            updateClothModel();

		Gui::DropdownList clothModelDropdown;
		clothModelDropdown.push_back({ ClothModel::FiniteElementsMethod, "Finite Elements Method" });
		clothModelDropdown.push_back({ ClothModel::ParticleSpringModel, "Particle Spring Model" });
		auto dirtyClothModel = pGui->addDropdown("Model", clothModelDropdown, reinterpret_cast<uint32_t&>(mClothModelType), true);

		dirtyClothModel |= pGui->addInt2Slider("Cloth Patch Tessellation", mClothPatchTessellation, 4, 128);
		dirtyClothModel |= pGui->addFloat2Slider("Cloth Patch Size", mClothPatchSize, 0.2f, 10.0f);

		if (dirtyClothModel)
			updateClothModel();

        if (pGui->addFloat3Slider("Wind Direction", mWindDirection, -1.0f, 1.0f))
        {
            if (length(mWindDirection) != 0.0f)
                mWindDirection = normalize(mWindDirection);
            showDirectionHelper(mWindDirection);
        }
        pGui->addFloatSlider("Wind Strength m/s", mWindStrength, 0.0f, 20.0f);
        pGui->addFloatSlider("Wind Drift Strength", mWindDriftStrength, 0.0f, 1.0f);

		pGui->endGroup();
	}

	mpClothModel->onGuiRender(this, pSample);
}

void ClothSample::onLoad(SampleCallbacks* pSample, RenderContext* pRenderContext)
{
    srand(uint32_t(time(0)));

	mpDbgUnitSphere = Model::createFromFile(ModelUnitSphere, Model::LoadFlags::None);
	mpDbgUnitCylinder = Model::createFromFile(ModelUnitCylinder, Model::LoadFlags::None);
	mpDbgUnitCone = Model::createFromFile(ModelUnitCone, Model::LoadFlags::None);

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

    mPhysicsObjects.emplace_back(new PhysicsSphere(vec3(0.0f, -2.0f, 4.0f), 1.0f));

	ClothModel::init();
	updateClothModel();

    Sampler::Desc samplerDesc;
    samplerDesc.setFilterMode(Sampler::Filter::Linear, Sampler::Filter::Linear, Sampler::Filter::Linear);
    mpTriLinearSampler = Sampler::create(samplerDesc);

    mpSkybox = SkyBox::create(SkyBoxTextures[linearRand<int32>(0, arraysize(SkyBoxTextures) - 1)], true, mpTriLinearSampler);

	mpScene = Scene::create();
	mpSceneRenderer = SceneRenderer::create(mpScene);
}

void ClothSample::onFrameRender(SampleCallbacks* pSample, RenderContext* pRenderContext, const Fbo::SharedPtr& pTargetFbo)
{
    const vec4 clearColor(0.38f, 0.52f, 0.10f, 1);
    pRenderContext->clearFbo(pTargetFbo.get(), clearColor, 1.0f, 0, FboAttachmentType::All);

	const auto currTime = pSample->getCurrentTime();
	const auto deltaTime = currTime - mPrevTime;
	mPrevTime = currTime;

	mpScene->deleteAllModels();
	simulateWind(deltaTime);
	mpClothModel->simulate(this, deltaTime);
    mpSkybox->render(pRenderContext, mpCamera.get());

    const auto showVectorField = [&](BoundingBox box, vec3 step, vec3 direction) -> void
    {
        const auto min = box.getMinPos();
        const auto max = box.getMaxPos();

        for (auto x = min.x; x <= max.x; x += step.x)
            for (auto y = min.y; y <= max.y; y += step.y)
                for (auto z = min.z; z <= max.z; z += step.z)
					drawVector(direction, vec3(x, y, z));
    };

    mpDirLight->setWorldDirection(mLightDirection);

    getActiveCameraController().update();
	mpClothModel->render(this, pSample);

    for (const auto &physicsObject : mPhysicsObjects)
        physicsObject->render(this, pSample);

    //BoundingBox box;
    //box.center = vec3(0, 0, 0);
    //box.extent = vec3(40, 40, 40);
    //showVectorField(box, vec3(10, 10, 10), mLightDirection);

	if (currTime < mHelperEndTime)
	{
		const auto &cameraData = mpCamera->getData();
		const auto position = vec3(glm::inverse(cameraData.viewMat) * vec4(0.9, -0.9 / cameraData.aspectRatio, -1.2f, 1.0f));
		drawVector(mHelperDirection, position);
	}

	mpSceneRenderer->update(currTime);
	mpDirLight->setIntoProgramVars(mpModelVars.get(), mpModelVars["PerFrameCB"].get(), "gDirLight");

	pRenderContext->setGraphicsVars(mpModelVars);
	pRenderContext->setGraphicsState(mpModelState);

	mpModelState->pushFbo(pTargetFbo);
	mpSceneRenderer->renderScene(pRenderContext, mpCamera.get());
    if (mpObjectController != nullptr)
        mpObjectController->renderScene(this, pSample);
	mpModelState->popFbo();
}

bool ClothSample::onKeyEvent(SampleCallbacks* pSample, const KeyboardEvent& keyEvent)
{
    return getActiveCameraController().onKeyEvent(keyEvent);
}

bool ClothSample::onMouseEvent(SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
    auto closestHit = std::numeric_limits<float>::infinity();
    switch (mouseEvent.type)
    {
    case MouseEvent::Type::LeftButtonDown:
    {
        SelectionQuery selectionQuery;
        GetMouseRay(mouseEvent, mpCamera->getData(), selectionQuery.ray.origin, selectionQuery.ray.direction);

        for (const auto &physicsObject : mPhysicsObjects)
            physicsObject->testSelection(selectionQuery);
        mpClothModel->testSelection(selectionQuery);

        if (selectionQuery.closestObject != nullptr)
        {
            if (mpSelectedObject != nullptr)
                mpSelectedObject->loseSelection(this);

            mpSelectedObject = selectionQuery.closestObject;
            mpSelectedObject->makeSelection(this, selectionQuery);
            return true;
        }
    } break;
    case MouseEvent::Type::LeftButtonUp:
    {
        if (mpSelectedObject != nullptr &&  mpSelectedObject->loseSelectionOnMouseUp())
        {
            mpSelectedObject->loseSelection(this);
            mpSelectedObject = nullptr;
            return true;
        }
    } break;
    }

    if (mpSelectedObject != nullptr && mpSelectedObject->onMouseEvent(this, pSample, mouseEvent))
        return true;

    return getActiveCameraController().onMouseEvent(mouseEvent);
}

void ClothSample::onResizeSwapChain(SampleCallbacks* pSample, uint32_t width, uint32_t height)
{
    mpCamera->setFocalLength(21.0f);
    float32_t aspectRatio = float32_t(width) / float32_t(height);
    mpCamera->setAspectRatio(aspectRatio);
}

vec3 ClothSample::getWindVelocity() const
{
    return mWindVelocity;
}

void ClothSample::simulateWind(float deltaTime)
{
    const auto &CurrentWindDrift = [&]() -> vec3 {
        return mWindDriftDirection * (mWindDriftStrength * sinf(pi<float>() * mWindDriftTimeNow / mWindDriftTimeTotal));
    };

    auto windDrift = vec3(0.0f);
    if (mWindDriftTimeNow < deltaTime)
    {
        windDrift = CurrentWindDrift();

        mWindDriftDirection = glm::sphericalRand(1.0f);
        mWindDriftTimeTotal = glm::linearRand(1.0f, 5.0f);
        mWindDriftTimeNow = mWindDriftTimeTotal + mWindDriftTimeNow;
    }

    mWindDriftTimeNow -= deltaTime;
    windDrift += CurrentWindDrift();
    mWindVelocity = (mWindDirection + windDrift) * mWindStrength;
}

void ClothSample::showDirectionHelper(const vec3 &vector, float time)
{
	mHelperDirection = normalize(vector);
	mHelperEndTime = mPrevTime + time;
}
void ClothSample::drawVector(const vec3 &vector, const vec3 &position)
{
	const auto yawPitchRoll = getYawPitchRoll(vector);

	auto unitSphere = Scene::ModelInstance::create(mpDbgUnitSphere, position, vec3(0.0f), vec3(0.0001f));
	unitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
	mpScene->addModelInstance(unitSphere);

	auto unitCylinder = Scene::ModelInstance::create(mpDbgUnitCylinder, position, yawPitchRoll, 0.01f * vec3(1.0f, length(vector), 1.0f));
	unitCylinder->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
	mpScene->addModelInstance(unitCylinder);

	auto unitCone = Scene::ModelInstance::create(mpDbgUnitCone, position + 0.1f * vector, yawPitchRoll, vec3(0.01f));
	unitCone->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
	mpScene->addModelInstance(unitCone);
}

bool ClothSample::addFloatSlider(const char label[], float& var, float minVal, float maxVal, bool sameLine, const char* displayFormat)
{
    if (sameLine) ImGui::SameLine();
    bool b = ImGui::SliderFloat(label, &var, minVal, maxVal, displayFormat);
    var = glm::clamp(var, minVal, maxVal);
    return b;
}

void GetMouseRay(const MouseEvent& mouseEvent, const CameraData &cameraData, vec3 &rayOrigin, vec3 &rayDirection)
{
    rayOrigin = cameraData.posW;
    rayDirection = cameraData.cameraW;
    rayDirection += cameraData.cameraU * (mouseEvent.pos.x - 0.5f) * 2.0f;
    rayDirection += cameraData.cameraV * (0.5f - mouseEvent.pos.y) * 2.0f;
    rayDirection = normalize(rayDirection);
}

int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nShowCmd)
{
    ClothSample::UniquePtr pRenderer = std::make_unique<ClothSample>();
    SampleConfig config;
    config.windowDesc.title = "Cloth Simulator";
    config.windowDesc.resizableWindow = true;
    config.deviceDesc.enableVsync = true;
    Sample::run(config, pRenderer);
    return 0;
}
