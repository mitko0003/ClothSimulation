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
#pragma once
#include "Falcor.h"
#include "ClothModel.h"

using namespace Falcor;

class ObjectController : std::enable_shared_from_this<ObjectController>
{
public:
    using SharedPtr = std::shared_ptr<ObjectController>;
    using SharedConstPtr = std::shared_ptr<const ObjectController>;

    ObjectController(ClothSample*, const glm::vec3 &translation, const glm::vec3 &scale);
    void renderScene(ClothSample*, SampleCallbacks*);
    void onGuiRender(ClothSample*, SampleCallbacks*);
    bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&);

    glm::vec3 getPosition() const;
    glm::vec3 getScale() const;

private:
    void activateGizmo(Gizmo::Type);

    Scene::SharedPtr mpScene;
    Picking::UniquePtr mpPicker;
    SceneRenderer::SharedPtr mpSceneRenderer;
    Gizmo::Gizmos mGizmos;
    Gizmo::Type mActiveGizmoType = Gizmo::Type::Translate;
    Scene::ModelInstance::SharedPtr mpProxyModelInstance;
    bool mGizmoBeingDragged = false;
    bool mUniformScale = true;
};

class PhysicsObject : public SelectableObject, std::enable_shared_from_this<PhysicsObject>
{
public:
    using SharedPtr = std::shared_ptr<PhysicsObject>;
    using SharedConstPtr = std::shared_ptr<const PhysicsObject>;

    virtual vec3 getPosition() const = 0;
    virtual vec3 getScale() const = 0;

    virtual float intersect(const Ray &ray) = 0;
    virtual vec3 intersectionResponse(vec3 V[3]) = 0;
    virtual void render(ClothSample*, SampleCallbacks*) = 0;
    virtual void applyController(const ObjectController*) = 0;

    // Physics object implementing selectable interface is ugly but we still do not have an entity component system to solve this issue.
    void testSelection(SelectionQuery&) override;
    void makeSelection(ClothSample*, SelectionQuery&) override;
    void loseSelection(ClothSample*) override;
    bool loseSelectionOnMouseUp() override { return false; }
    bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) override;

private:
    ObjectController::SharedPtr mpController;
};

class ClothSample : public Renderer
{
public:
    void onLoad(SampleCallbacks* pSample, RenderContext* pRenderContext) override;
    void onFrameRender(SampleCallbacks* pSample, RenderContext* pRenderContext, const Fbo::SharedPtr& pTargetFbo) override;
	void onShutdown(SampleCallbacks* pSample) override {}
    void onResizeSwapChain(SampleCallbacks* pSample, uint32_t width, uint32_t height) override;
    bool onKeyEvent(SampleCallbacks* pSample, const KeyboardEvent& keyEvent) override;
    bool onMouseEvent(SampleCallbacks* pSample, const MouseEvent& mouseEvent) override;
	void onDataReload(SampleCallbacks* pSample) override {}
    void onGuiRender(SampleCallbacks* pSample, Gui* pGui) override;
	void showDirectionHelper(const vec3&, float time = 2.0f);
	void drawVector(const vec3 &vector, const vec3 &position);
    vec3 getWindVelocity() const;

    // Gui::addFloatSlider ignores displayFormat. Use this function until a Falcor version with fix comes.
    static bool addFloatSlider(const char label[], float& var, float minVal, float maxVal, bool sameLine, const char* displayFormat = "%.3f");

	Scene::SharedPtr mpScene;
	SceneRenderer::SharedPtr mpSceneRenderer;

	Camera::SharedPtr mpCamera;
	DirectionalLight::SharedPtr mpDirLight;

	GraphicsProgram::SharedPtr mpModelProgram;
	GraphicsVars::SharedPtr mpModelVars;
	GraphicsState::SharedPtr mpModelState;

	Model::SharedPtr mpDbgUnitSphere;
	Model::SharedPtr mpDbgUnitCylinder;
	Model::SharedPtr mpDbgUnitCone;

    std::vector<PhysicsObject::SharedPtr> mPhysicsObjects;
    ObjectController::SharedPtr mpObjectController;

private:
    static const char *ModelPS;
	static const char *ModelUnitSphere;
	static const char *ModelUnitCylinder;
	static const char *ModelUnitCone;

    static const char *SkyBoxTextures[];

    SkyBox::SharedPtr mpSkybox;
    Sampler::SharedPtr mpTriLinearSampler;

    enum
    {
        ModelViewCamera,
        FirstPersonCamera,
        SixDoFCamera
    } mCameraType;

    ModelViewCameraController mModelViewCameraController;
    FirstPersonCameraController mFirstPersonCameraController;
    SixDoFCameraController mSixDoFCameraController;

    CameraController& getActiveCameraController();
	void updateClothModel();

    vec3 mWindVelocity = vec3(0.0f);
    vec3 mWindDirection = vec3(1.0f, 0.0f, 0.0f);
    vec3 mWindDriftDirection = vec3(0.0f);
    float mWindDriftTimeTotal = 1.0f;
    float mWindDriftTimeNow = 0.0f;
    float mWindDriftStrength = 0.5f;
    float mWindStrength = 0.0f;
    void simulateWind(float deltaTime);

	vec3 mHelperDirection;
    vec3 mLightDirection;

    float mPrevTime;
	float mHelperEndTime;

	vec2 mClothPatchSize = vec2(5.0f);
	ivec2 mClothPatchTessellation = ivec2(12, 12);
	ClothModel::EType mClothModelType = ClothModel::ParticleSpringModel;

    SelectableObject *mpSelectedObject;

    ClothModel::SharedPtr mpClothModel;
};

void GetMouseRay(const MouseEvent& mouseEvent, const CameraData &cameraData, vec3 &rayOrigin, vec3 &rayDirection);
