#pragma once

#include "Falcor.h"
#include "Utils.h"
#include "External/HalfEdge/trimesh.h"

using namespace Falcor;
struct SelectionQuery;
class ClothSample;

struct SelectableObject
{
    virtual ~SelectableObject() {}

    virtual void testSelection(SelectionQuery&) = 0;
    virtual void makeSelection(ClothSample*, SelectionQuery&) = 0;
    virtual void loseSelection(ClothSample*) = 0;

    virtual bool loseSelectionOnMouseUp() { return true; }
    virtual bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) = 0;
};

struct SelectionQuery
{
    Ray ray;
    float closestHit = std::numeric_limits<float>::infinity();
    SelectableObject *closestObject = nullptr;
    uint8 cache[256]; // owned by current closest object
};

struct ClothModel : NonCopyable, public SelectableObject, public std::enable_shared_from_this<ClothModel>
{
	using SharedPtr = std::shared_ptr<ClothModel>;
	using SharedConstPtr = std::shared_ptr<const ClothModel>;

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

	Buffer::SharedPtr mpVBPositions;
	Buffer::SharedPtr mpVBNormals;
	Buffer::SharedPtr mpVBColors;
	Buffer::SharedPtr mpIndexBuffer;

	Vao::SharedPtr mpVao;

	GraphicsState::SharedPtr mpClothState;
	RasterizerState::SharedPtr mpRasterizeNormal;
	RasterizerState::SharedPtr mpRasterizeWireframe;

	virtual EType getType() const = 0;
	virtual vec3 getVertexPosition(int32 index) const = 0;
    virtual vec3 getVertexNormal(int32 index) const;
    void getTriangleVertices(const trimesh::triangle_t &triangle, vec3 vertices[3]) const;
    vec3 computeVertexNormal(int32 index) const;

	virtual void init(const vec2 &size, const ivec2 &tessellation) = 0;
	virtual void init(const Model *model) = 0;
	virtual void simulate(ClothSample*, float deltaTime) = 0;
	virtual void render(ClothSample*, SampleCallbacks*);
	virtual void onGuiRender(ClothSample*, SampleCallbacks*);

    virtual void testSelection(SelectionQuery&) override {};
    virtual void makeSelection(ClothSample*, SelectionQuery&) override {};
    virtual void loseSelection(ClothSample*) override {};

    virtual bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) { return false; }

	static void init();
	static int32 createRectMesh(const vec2 &size, const ivec2 &tessellation, std::vector<float3> &position, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh);

	static SharedPtr createClothModel(EType type, const ClothModel *parent);

protected:
	void sharedInit();

	struct
	{
        float windDragCoefficient = 0.1f;
        float windLiftCoefficient = 0.1f;
		bool bShowWindEffect = false;
		bool bShowNormals = false;
		bool bShowParticles = false;
		bool bShowWireframe = false;
		bool bShowSelfCollisions = false;
        vec3 baseColor = vec3(0.0f, 1.0f, 0.0f);
	} mUserParams;

	template <EType type>
	static SharedPtr createClothModel(const ClothModel *parent);
};
