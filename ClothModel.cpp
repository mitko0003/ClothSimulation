#include "ClothModel.h"

const char *ClothModel::ModelUnitSphere = "Models/UnitSphere.fbx";
const char *ClothModel::ModelUnitCylinder = "Models/UnitCylinder.fbx";
const char *ClothModel::ModelUnitCone = "Models/UnitCone.fbx";

int32 ClothModel::createRectMesh(int32_t sizeX, int32_t sizeY, std::vector<float3> &positions, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh)
{
	const auto vertexCount = sizeX * sizeY;

	for (int32_t y = 0; y < sizeY; ++y)
	{
		for (int32_t x = 0; x < sizeX; ++x)
		{
			positions.emplace_back(float32(x), float32(y), 0.0f);
		}
	}

	const auto EmitTriangle = [&](const auto &v0, const auto &v1, const auto &v2) -> void
	{
		trimesh::triangle_t triangle;
		triangle.v[0] = v0;
		triangle.v[1] = v1;
		triangle.v[2] = v2;
		triangles.push_back(triangle);
	};

	for (int32_t y = 0; y < sizeY - 1; ++y)
	{
		for (int32_t x = 0; x < sizeX - 1; ++x)
		{
			EmitTriangle(
				(y + 1) * sizeY + x,
				y * sizeY + (x + 1),
				y * sizeY + x
			);
			EmitTriangle(
				(y + 1) * sizeY + (x + 1),
				y * sizeY + (x + 1),
				(y + 1) * sizeY + x
			);
		}
	}

	std::vector<trimesh::edge_t> edges;
	trimesh::unordered_edges_from_triangles(triangles.size(), &triangles[0], edges);

	mesh.build(vertexCount, triangles.size(), &triangles[0], edges.size(), &edges[0]);

	return vertexCount;
}