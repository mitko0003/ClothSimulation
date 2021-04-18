#include <algorithm>
#include "TriangleIntersection.h"
#include "ClothSample.h"
#include "ClothModel.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#pragma optimize("", off)
const char *ClothModel::ClothVS = "Cloth.vs.hlsl";
const char *ClothModel::ClothPS = "Cloth.ps.hlsl";

GraphicsProgram::SharedPtr ClothModel::spClothProgram = nullptr;
GraphicsVars::SharedPtr ClothModel::spClothVars = nullptr;
VertexLayout::SharedPtr ClothModel::spVertexLayout = nullptr;
std::vector<ClothModel::TTextureSet> ClothModel::sTextureSets;

int32 ClothModel::createRectMesh(const vec2 &size, const ivec2 &tessellation, std::vector<float3> &positions, std::vector<float2> &texCoords, std::vector<trimesh::triangle_t> &triangles, trimesh::trimesh_t &mesh)
{
	const auto vertexCount = tessellation.x * tessellation.y;
	const auto gridStep = size / vec2(tessellation - 1);

	for (int32 y = 0; y < tessellation.y; ++y)
	{
		for (int32 x = 0; x < tessellation.x; ++x)
		{
			positions.emplace_back(float32(x) * gridStep.x - size.x * 0.5f, size.y * 0.5f - float32(y) * gridStep.y, 0.0f);
            texCoords.emplace_back(float2(x, y) / float2(tessellation - 1));
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

	for (int32 y = 0; y < tessellation.y - 1; ++y)
	{
		for (int32 x = 0; x < tessellation.x - 1; ++x)
		{
			EmitTriangle(
				(y + 1) * tessellation.x + x,
				y * tessellation.x + (x + 1),
				y * tessellation.x + x
			);
			EmitTriangle(
				(y + 1) * tessellation.x + (x + 1),
				y * tessellation.x + (x + 1),
				(y + 1) * tessellation.x + x
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
    VertexBufferLayout::SharedPtr pTexCoordsLayout = VertexBufferLayout::create();
    pTexCoordsLayout->addElement("TEXCOORD", 0, ResourceFormat::RG16Unorm, 1, 3);
    spVertexLayout->addBufferLayout(3, pTexCoordsLayout);

    for (auto& file : fs::directory_iterator("Data/Materials/"))
    {
        const auto &ReadTexture = [=](const char *name, bool generateMipLevels, bool loadAsSrgb) -> Texture::SharedPtr
        {
            auto path = file.path().string() + "/" + name;
            return createTextureFromFile(path, generateMipLevels, loadAsSrgb, Texture::BindFlags::ShaderResource);
        };
        std::string filename = file.path().filename().string();

        TTextureSet textureSet;
        textureSet.Name = file.path().filename().string();
        textureSet.BaseColorMap = ReadTexture("BaseColor.png", true, true);
        textureSet.RoughnessMap = ReadTexture("Roughness.png", true, false);
        textureSet.NormalMap = ReadTexture("Normal.png", true, false);
        sTextureSets.emplace_back(std::move(textureSet));
    }

	Program::DefineList defineList;

	GraphicsProgram::Desc clothProgramDesc;
	clothProgramDesc.addShaderLibrary(ClothVS).vsEntry("main");
	clothProgramDesc.addShaderLibrary(ClothPS).psEntry("main");
	spClothProgram = GraphicsProgram::create(clothProgramDesc, defineList);
	spClothVars = GraphicsVars::create(spClothProgram->getReflector());
}

void ClothModel::sharedInit(const std::vector<float2> &texCoords)
{
	std::vector<uint16> vIndices;
	for (const auto &triangle : mTriangles)
	{
		vIndices.emplace_back(uint16(triangle.v[0]));
		vIndices.emplace_back(uint16(triangle.v[1]));
		vIndices.emplace_back(uint16(triangle.v[2]));
	}

    std::vector<uint16> unormTexCoords;
    for (const auto &uv : texCoords)
    {
        unormTexCoords.push_back(uint16_t(uv.x * std::numeric_limits<uint16_t>::max()));
        unormTexCoords.push_back(uint16_t(uv.y * std::numeric_limits<uint16_t>::max()));
    }

	const auto vertexCount = mTriangles.size() * 3;

	mpIndexBuffer = Buffer::create(vIndices.size() * sizeof(uint16_t), ResourceBindFlags::Index, Buffer::CpuAccess::None, vIndices.data());
	mpVBPositions = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBNormals = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBColors = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBTexCoords = Buffer::create(vertexCount * sizeof(uint16_t[2]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, unormTexCoords.data());
    
	Vao::BufferVec buffers{ mpVBPositions, mpVBNormals, mpVBColors, mpVBTexCoords };
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

void ClothModel::getTriangleVertices(const trimesh::triangle_t &triangle, vec3 vertices[3]) const
{
    vertices[0] = getVertexPosition(triangle[0]);
    vertices[1] = getVertexPosition(triangle[1]);
    vertices[2] = getVertexPosition(triangle[2]);
}

vec3 ClothModel::getVertexNormal(int32 index) const
{
    return computeVertexNormal(index);
}

vec3 ClothModel::computeVertexNormal(int32 index) const
{
    auto faces = mMesh.vertex_face_neighbors(index);
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

    return normalize(normal);
}

void ClothModel::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	const auto currTime = pSample->getCurrentTime();
	auto *pRenderContext = pSample->getRenderContext();
	const auto &pTargetFbo = pSample->getCurrentFbo();

	const int32 vertexCount = int32(mMesh.get_num_vertices());

	float32_t *pPositions = reinterpret_cast<float32_t*>(mpVBPositions->map(Buffer::MapType::WriteDiscard));
	float32_t *pNormals = reinterpret_cast<float32_t*>(mpVBNormals->map(Buffer::MapType::WriteDiscard));
	for (int32 vertex = 0; vertex < vertexCount; ++vertex)
	{
		const auto position = getVertexPosition(vertex);
		std::copy_n(position.data.data, 3, pPositions + vertex * 3);

        const auto normal = getVertexNormal(vertex);
		std::copy_n(normal.data.data, 3, pNormals + vertex * 3);

		if (mUserParams.bShowNormals)
		{
			pClothSample->drawVector(normal, position);
		}
		else if (mUserParams.bShowParticles)
		{
			auto unitSphere = Scene::ModelInstance::create(pClothSample->mpDbgUnitSphere, position, vec3(0.0f), vec3(0.0005f));
			unitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
			pClothSample->mpScene->addModelInstance(unitSphere);
		}
	}
	mpVBPositions->unmap();
	mpVBNormals->unmap();

	auto *pColors = reinterpret_cast<vec3*>(mpVBColors->map(Buffer::MapType::WriteDiscard));
	std::fill_n(pColors, vertexCount, mUserParams.baseColor);

	if (mUserParams.bShowSelfCollisions)
	{
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
            getTriangleVertices(triangle0, t0_vertices);

			for (int32 t1 = t0 + 1; t1 < mTriangles.size(); ++t1)
			{
				vec3 t1_vertices[3];
				const auto &triangle1 = mTriangles[t1];
                getTriangleVertices(triangle1, t1_vertices);

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
	}
	mpVBColors->unmap();

	if (mUserParams.bShowWireframe)
		mpClothState->setRasterizerState(mpRasterizeWireframe);
	else
		mpClothState->setRasterizerState(mpRasterizeNormal);

	pClothSample->mpCamera->setIntoConstantBuffer(spClothVars->getConstantBuffer("InternalPerFrameCB").get(), "gCamera");
	pClothSample->mpDirLight->setIntoProgramVars(spClothVars.get(), spClothVars["PerFrameCB"].get(), "gDirLight");

    if (!sTextureSets.empty())
    {
        const auto &textureSet = sTextureSets[mUserParams.textureSet];
        spClothVars->setTexture("gBaseColorMap", textureSet.BaseColorMap);
        spClothVars->setTexture("gRoughnessMap", textureSet.RoughnessMap);
        spClothVars->setTexture("gNormalMap", textureSet.NormalMap);
    }

	pRenderContext->setGraphicsVars(spClothVars);
	pRenderContext->setGraphicsState(mpClothState);

	mpClothState->pushFbo(pTargetFbo);
	pRenderContext->drawIndexed(uint32(mTriangles.size()) * 3u, 0, 0);
	mpClothState->popFbo();
}

void ClothModel::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

    pGui->addFloatSlider("Drag", mUserParams.windDragCoefficient, 0.0f, 1.0f);
    pGui->addFloatSlider("Lift", mUserParams.windLiftCoefficient, 0.0f, 1.0f);

	pGui->addCheckBox("Show Particles", mUserParams.bShowParticles);
	pGui->addCheckBox("Show Wireframe", mUserParams.bShowWireframe);
    pGui->addCheckBox("Show Surface Normals", mUserParams.bShowNormals);
    pGui->addCheckBox("Show Wind Effect", mUserParams.bShowWindEffect);
	pGui->addCheckBox("Show Self Collisions", mUserParams.bShowSelfCollisions);

    if (!sTextureSets.empty())
    {
        Gui::DropdownList textureSetDropdown;
        for (uint32 textureSet = 0u; textureSet < sTextureSets.size(); ++textureSet)
            textureSetDropdown.push_back({ textureSet, sTextureSets[textureSet].Name });
        pGui->addDropdown("Texture Set", textureSetDropdown, reinterpret_cast<uint32&>(mUserParams.textureSet));
    }

    pGui->addRgbColor("Base Color", mUserParams.baseColor);
}

auto ClothModel::createClothModel(EType type, const ClothModel *parent) -> SharedPtr
{
	const auto *parentToPass = parent && parent->getType() == type ? parent : nullptr;
	SharedPtr result;
	switch (type)
	{
	case ParticleSpringModel: result = createClothModel<ParticleSpringModel>(parentToPass); break;
	case FiniteElementsMethod: result = createClothModel<FiniteElementsMethod>(parentToPass); break;
	default:
		should_not_get_here();
		return nullptr;
	}

	if (parent != nullptr)
	{
		result->mUserParams = parent->mUserParams;
	}

	return result;
}
