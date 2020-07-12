#pragma once

#include "Falcor.h"
#include "Utils.h"
#include "External/HalfEdge/trimesh.h"

using namespace Falcor;

class ClothSample;

struct ClothModel : NonCopyable, std::enable_shared_from_this<ClothModel>
{
	static const char *ModelUnitSphere;
	static const char *ModelUnitCylinder;
	static const char *ModelUnitCone;

	std::string mName;

	virtual void init(int32_t width, int32_t height) = 0;
	virtual void init(const Model *model) = 0;
	virtual void simulate(float deltaTime) = 0;
	virtual void render(ClothSample*, SampleCallbacks*) = 0;

	virtual void onGuiRender(ClothSample*, SampleCallbacks*) = 0;
	virtual bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) = 0;

	static void init();
	static int32 createRectMesh(int32_t width, int32_t height, std::vector<float3> &position, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh);
};