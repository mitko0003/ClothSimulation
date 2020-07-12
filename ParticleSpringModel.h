#pragma once

#include "ClothModel.h"
#include "External/HalfEdge/trimesh.h"

struct ClothParticle
{
	float3 mPosition;
	float3 mPrevPosition;
	float3 mAcceleration;
	bool mbStationary;

	Scene::ModelInstance::SharedPtr pDbgUnitSphere;
	Scene::ModelInstance::SharedPtr pDbgUnitCylinder;
};

struct ClothPatch : public ClothModel
{
	static const char *ClothVS;
	static const char *ClothPS;

	ClothParticle *mpParticles;
	int32_t mWidth;
	int32_t mHeight;
	float32 mMass;

	std::vector<trimesh::triangle_t> triangles;
	trimesh::trimesh_t mesh;

	Buffer::SharedPtr mpVBPositions;
	Buffer::SharedPtr mpVBNormals;
	Buffer::SharedPtr mpVBColors;
	Buffer::SharedPtr mpIndexBuffer;

	ClothParticle *mHeldParticle;
	float32 mAirTemperature;

	bool mbShowParticles;
	bool mbShowWireframe;

	GraphicsProgram::SharedPtr mpClothProgram;
	GraphicsVars::SharedPtr mpClothVars;
	GraphicsState::SharedPtr mpClothState;

	RasterizerState::SharedPtr mpRasterizeNormal;
	RasterizerState::SharedPtr mpRasterizeWireframe;

	Vao::SharedPtr mpVao;

	Scene::SharedPtr mpDbgScene;
	SceneRenderer::SharedPtr mpDbgSceneRenderer;

	ClothPatch() :
		mHeldParticle(nullptr),
		mpParticles(nullptr),
		mWidth(0),
		mHeight(0),
		mMass(7.0f)
	{
		mName = "Varlet Cloth Patch";
	}

	~ClothPatch()
	{
		delete[] mpParticles;
	}

	ClothParticle* getParticle(int32_t x, int32_t y);

	void init(int32_t width, int32_t height) final override;
	void init(const Model *model) final override;
	void simulate(float deltaTime) final override;
	void render(ClothSample*, SampleCallbacks*) final override;

	void onGuiRender(ClothSample*, SampleCallbacks*) final override;
	bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) final override;
};