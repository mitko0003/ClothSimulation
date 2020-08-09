#include <algorithm>
#include "TriangleIntersection.h"
#include "ClothSample.h"
#include "ClothModel.h"

const char *ClothModel::ClothVS = "Cloth.vs.hlsl";
const char *ClothModel::ClothPS = "Cloth.ps.hlsl";

GraphicsProgram::SharedPtr ClothModel::spClothProgram = nullptr;
GraphicsVars::SharedPtr ClothModel::spClothVars = nullptr;
VertexLayout::SharedPtr ClothModel::spVertexLayout = nullptr;

int32 ClothModel::createRectMesh(int32_t sizeX, int32_t sizeY, std::vector<float3> &positions, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh)
{
	const auto vertexCount = sizeX * sizeY;

	for (int32 y = 0; y < sizeY; ++y)
	{
		for (int32 x = 0; x < sizeX; ++x)
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

	for (int32 y = 0; y < sizeY - 1; ++y)
	{
		for (int32 x = 0; x < sizeX - 1; ++x)
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
	trimesh::unordered_edges_from_triangles(int(triangles.size()), &triangles[0], edges);

	mesh.build(vertexCount, int(triangles.size()), &triangles[0], int(edges.size()), &edges[0]);

	return vertexCount;
}

void ClothModel::init()
{
	spVertexLayout = VertexLayout::create();
	VertexBufferLayout::SharedPtr pPositionsLayout = VertexBufferLayout::create();
	pPositionsLayout->addElement("POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
	spVertexLayout->addBufferLayout(0, pPositionsLayout);
	VertexBufferLayout::SharedPtr pNormalsLayout = VertexBufferLayout::create();
	pNormalsLayout->addElement("NORMAL", 0, ResourceFormat::RGB32Float, 1, 1);
	spVertexLayout->addBufferLayout(1, pNormalsLayout);
	VertexBufferLayout::SharedPtr pColorsLayout = VertexBufferLayout::create();
	pColorsLayout->addElement("COLOR", 0, ResourceFormat::RGB32Float, 1, 2);
	spVertexLayout->addBufferLayout(2, pColorsLayout);

	Program::DefineList defineList;

	GraphicsProgram::Desc clothProgramDesc;
	clothProgramDesc.addShaderLibrary(ClothVS).vsEntry("main");
	clothProgramDesc.addShaderLibrary(ClothPS).psEntry("main");
	spClothProgram = GraphicsProgram::create(clothProgramDesc, defineList);
	spClothVars = GraphicsVars::create(spClothProgram->getReflector());
}

void ClothModel::sharedInit()
{
	std::vector<uint16> vIndices;
	for (const auto &triangle : mTriangles)
	{
		vIndices.emplace_back(uint16(triangle.v[0]));
		vIndices.emplace_back(uint16(triangle.v[1]));
		vIndices.emplace_back(uint16(triangle.v[2]));
	}

	const auto vertexCount = mTriangles.size() * 3;

	mpIndexBuffer = Buffer::create(vIndices.size() * sizeof(uint16_t), ResourceBindFlags::Index, Buffer::CpuAccess::None, vIndices.data());
	mpVBPositions = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBNormals = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBColors = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);

	Vao::BufferVec buffers{ mpVBPositions, mpVBNormals, mpVBColors };
	mpVao = Vao::create(Vao::Topology::TriangleList, spVertexLayout, buffers, mpIndexBuffer, ResourceFormat::R16Uint);

	mpClothState = GraphicsState::create();

	RasterizerState::Desc rasterizerDesc;
	rasterizerDesc.setCullMode(RasterizerState::CullMode::None);
	rasterizerDesc.setFillMode(RasterizerState::FillMode::Solid);
	mpRasterizeNormal = RasterizerState::create(rasterizerDesc);
	rasterizerDesc.setFillMode(RasterizerState::FillMode::Wireframe);
	mpRasterizeWireframe = RasterizerState::create(rasterizerDesc);

	mpClothState->setProgram(spClothProgram);
	mpClothState->setVao(mpVao);
}

void ClothModel::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	const auto currTime = pSample->getCurrentTime();
	auto *pRenderContext = pSample->getRenderContext();
	const auto &pTargetFbo = pSample->getCurrentFbo();

	Scene::SharedPtr dbgScene =
		mbShowParticles ? Scene::create() : nullptr;
	SceneRenderer::SharedPtr dbgSceneRenderer =
		mbShowParticles ? SceneRenderer::create(dbgScene) : nullptr;

	const int32 vertexCount = int32(mMesh.get_num_vertices());

	float32_t *pPositions = reinterpret_cast<float32_t*>(mpVBPositions->map(Buffer::MapType::WriteDiscard));
	float32_t *pNormals = reinterpret_cast<float32_t*>(mpVBNormals->map(Buffer::MapType::WriteDiscard));
	for (int32 vertex = 0; vertex < vertexCount; ++vertex)
	{
		const auto position = getVertexPosition(vertex);
		std::copy_n(position.data.data, 3, pPositions + vertex * 3);

		auto faces = mMesh.vertex_face_neighbors(vertex);
		float3 normal(0.0f);

		const auto &GetPlaneNormal = [](const float3 &A, const float3 &B, const float3 &C) -> float3 {
			const auto result = normalize(cross(C - B, A - B));
			return any(isnan(result)) ? float3(0.0f) : result;
		};

		for (const auto &face : faces)
		{
			const auto &triangle = mTriangles[face];
			normal += GetPlaneNormal(
				getVertexPosition(triangle[0]),
				getVertexPosition(triangle[1]),
				getVertexPosition(triangle[2])
			);
		}

		normal = normalize(normal);
		std::copy_n(normal.data.data, 3, pNormals + vertex * 3);

		if (mbShowParticles)
		{
			auto unitSphere = Scene::ModelInstance::create(pClothSample->mpDbgUnitSphere, position, vec3(0.0f), vec3(0.1f));
			unitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
			dbgScene->addModelInstance(unitSphere);
		}
	}
	mpVBPositions->unmap();
	mpVBNormals->unmap();

	const auto &GetTriangleVertices = [&](const trimesh::triangle_t &triangle, vec3 vertices[3]) -> void
	{
		vertices[0] = getVertexPosition(triangle[0]);
		vertices[1] = getVertexPosition(triangle[1]);
		vertices[2] = getVertexPosition(triangle[2]);
	};

	auto *pColors = reinterpret_cast<vec3*>(mpVBColors->map(Buffer::MapType::WriteDiscard));
	std::fill_n(pColors, vertexCount, vec3(0.0f, 1.0f, 0.0f));

	const auto &HandleCollision = [&](const trimesh::triangle_t &triangle) -> void
	{
		pColors[triangle[0]] = vec3(1.0f, 0.0f, 0.0f);
		pColors[triangle[1]] = vec3(1.0f, 0.0f, 0.0f);
		pColors[triangle[2]] = vec3(1.0f, 0.0f, 0.0f);
		// getParticle(x, y)->mPosition = getParticle(x, y)->mPrevPosition;
	};

	for (int32 t0 = 0; t0 < mTriangles.size(); ++t0)
	{
		vec3 t0_vertices[3];
		const auto &triangle0 = mTriangles[t0];
		GetTriangleVertices(triangle0, t0_vertices);

		for (int32 t1 = t0 + 1; t1 < mTriangles.size(); ++t1)
		{
			vec3 t1_vertices[3];
			const auto &triangle1 = mTriangles[t1];
			GetTriangleVertices(triangle1, t1_vertices);

			vec3 response; // TODO(mitko): Remove extra "responce" code from triangle intersect implementation
			bool collision = intersectTriangles(t0_vertices, t1_vertices, response);
			if (!collision) continue;

			for (const auto vertex : triangle0.v)
			{
				const auto faces = mMesh.vertex_face_neighbors(vertex);
				for (const auto face : faces)
				{
					collision &= triangle1[0] != vertex;
					collision &= triangle1[1] != vertex;
					collision &= triangle1[2] != vertex;
				}
			}

			if (!collision) continue;
			HandleCollision(triangle0);
			HandleCollision(triangle1);
		}
	}
	mpVBColors->unmap();

	if (mbShowWireframe)
		mpClothState->setRasterizerState(mpRasterizeWireframe);
	else
		mpClothState->setRasterizerState(mpRasterizeNormal);

	pClothSample->mpCamera->setIntoConstantBuffer(spClothVars->getConstantBuffer("InternalPerFrameCB").get(), "gCamera");
	pClothSample->mpDirLight->setIntoProgramVars(spClothVars.get(), spClothVars["PerFrameCB"].get(), "gDirLight");

	pRenderContext->setGraphicsVars(spClothVars);
	pRenderContext->setGraphicsState(mpClothState);

	mpClothState->pushFbo(pTargetFbo);
	pRenderContext->drawIndexed(uint32(mTriangles.size()) * 3u, 0, 0);
	mpClothState->popFbo();

	if (dbgSceneRenderer != nullptr)
	{
		dbgSceneRenderer->update(currTime);
		pClothSample->mpDirLight->setIntoProgramVars(
			pClothSample->mpModelVars.get(), pClothSample->mpModelVars["PerFrameCB"].get(), "gDirLight"
		);

		pRenderContext->setGraphicsVars(pClothSample->mpModelVars);
		pRenderContext->setGraphicsState(pClothSample->mpModelState);

		pClothSample->mpModelState->pushFbo(pTargetFbo);
		dbgSceneRenderer->renderScene(pRenderContext, pClothSample->mpCamera.get());
		pClothSample->mpModelState->popFbo();
	}
}

void ClothModel::onGuiRender(ClothSample*, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();
	pGui->addCheckBox("Show Particles", mbShowParticles);
	pGui->addCheckBox("Show Wireframe", mbShowWireframe);
}

ClothModel *ClothModel::createClothModel(EType type)
{
	switch (type)
	{
	case ParticleSpringModel: return createClothModel<ParticleSpringModel>();
	case FiniteElementsMethod: return createClothModel<FiniteElementsMethod>();
	default:
		should_not_get_here();
		return nullptr;
	}
}