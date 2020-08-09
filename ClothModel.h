#pragma once

#include "Falcor.h"
#include "Utils.h"
#include "External/HalfEdge/trimesh.h"

using namespace Falcor;

class ClothSample;

struct ClothModel : NonCopyable, std::enable_shared_from_this<ClothModel>
{
	enum EType
	{
		FiniteElementsMethod,
		ParticleSpringModel
	};

	static const char *ClothVS;
	static const char *ClothPS;

	static GraphicsProgram::SharedPtr spClothProgram;
	static GraphicsVars::SharedPtr spClothVars;
	static VertexLayout::SharedPtr spVertexLayout;

	std::string mName;

	std::vector<trimesh::triangle_t> mTriangles;
	trimesh::trimesh_t mMesh;

	float32 mMaterialDensity = 1.0f;

	Buffer::SharedPtr mpVBPositions;
	Buffer::SharedPtr mpVBNormals;
	Buffer::SharedPtr mpVBColors;
	Buffer::SharedPtr mpIndexBuffer;

	Vao::SharedPtr mpVao;

	GraphicsState::SharedPtr mpClothState;
	RasterizerState::SharedPtr mpRasterizeNormal;
	RasterizerState::SharedPtr mpRasterizeWireframe;

	virtual EType getType() const = 0;
	virtual glm::vec3 getVertexPosition(int32 index) const = 0;

	virtual void init(int32_t width, int32_t height) = 0;
	virtual void init(const Model *model) = 0;
	virtual void simulate(float deltaTime) = 0;
	virtual void render(ClothSample*, SampleCallbacks*);

	virtual void onGuiRender(ClothSample*, SampleCallbacks*);
	virtual bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) = 0;

	static void init();
	static int32 createRectMesh(int32_t width, int32_t height, std::vector<float3> &position, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh);

	static ClothModel *createClothModel(EType type);

protected:
	void sharedInit();

	bool mbShowParticles = false;
	bool mbShowWireframe = false;

	template <EType type>
	static ClothModel *createClothModel();
};