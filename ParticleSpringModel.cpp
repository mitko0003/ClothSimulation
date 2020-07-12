#include "ClothSample.h"
#include "ParticleSpringModel.h"
#include "TriangleIntersection.h"

using float1 = float32_t;
using float3 = f32vec3;

const char *ClothPatch::ClothVS = "Cloth.vs.hlsl";
const char *ClothPatch::ClothPS = "Cloth.ps.hlsl";

void ClothPatch::init(int32_t width, int32_t height)
{
	std::vector<float3> positions;
	const auto vertexCount = createRectMesh(width, height, positions, triangles, mesh);

	mWidth = width;
	mHeight = height;
	mpParticles = new ClothParticle[vertexCount];

	for (int32_t i = 0; i < vertexCount; ++i)
	{
		auto &particle = mpParticles[i];
		particle.mAcceleration = float3(0.0f, 0.0f, 0.0f);
		particle.mbStationary = false;
		particle.mPosition = positions[i];
		particle.mPrevPosition = particle.mPosition;
	}

	mpVBPositions = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBNormals = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);
	mpVBColors = Buffer::create(vertexCount * sizeof(float32_t[3]), ResourceBindFlags::Vertex, Buffer::CpuAccess::Write, nullptr);

	std::vector<uint16_t> vIndices;
	for (const auto &triangle : triangles)
	{
		vIndices.emplace_back(uint16_t(triangle.v[0]));
		vIndices.emplace_back(uint16_t(triangle.v[1]));
		vIndices.emplace_back(uint16_t(triangle.v[2]));
	}
	mpIndexBuffer = Buffer::create(vIndices.size() * sizeof(uint16_t), ResourceBindFlags::Index, Buffer::CpuAccess::None, vIndices.data());

	VertexLayout::SharedPtr pLayout = VertexLayout::create();
	VertexBufferLayout::SharedPtr pPositionsLayout = VertexBufferLayout::create();
	pPositionsLayout->addElement("POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
	pLayout->addBufferLayout(0, pPositionsLayout);
	VertexBufferLayout::SharedPtr pNormalsLayout = VertexBufferLayout::create();
	pNormalsLayout->addElement("NORMAL", 0, ResourceFormat::RGB32Float, 1, 1);
	pLayout->addBufferLayout(1, pNormalsLayout);
	VertexBufferLayout::SharedPtr pColorsLayout = VertexBufferLayout::create();
	pColorsLayout->addElement("COLOR", 0, ResourceFormat::RGB32Float, 1, 2);
	pLayout->addBufferLayout(2, pColorsLayout);

	Vao::BufferVec buffers{ mpVBPositions, mpVBNormals, mpVBColors };
	mpVao = Vao::create(Vao::Topology::TriangleList, pLayout, buffers, mpIndexBuffer, ResourceFormat::R16Uint);

	Program::DefineList defineList;

	GraphicsProgram::Desc clothProgramDesc;
	clothProgramDesc.addShaderLibrary(ClothVS).vsEntry("main");
	clothProgramDesc.addShaderLibrary(ClothPS).psEntry("main");
	mpClothProgram = GraphicsProgram::create(clothProgramDesc, defineList);
	mpClothVars = GraphicsVars::create(mpClothProgram->getReflector());
	mpClothState = GraphicsState::create();

	RasterizerState::Desc rasterizerDesc;
	rasterizerDesc.setCullMode(RasterizerState::CullMode::None);
	rasterizerDesc.setFillMode(RasterizerState::FillMode::Solid);
	mpRasterizeNormal = RasterizerState::create(rasterizerDesc);
	rasterizerDesc.setFillMode(RasterizerState::FillMode::Wireframe);
	mpRasterizeWireframe = RasterizerState::create(rasterizerDesc);

	mpClothState->setProgram(mpClothProgram);
	mpClothState->setVao(mpVao);

	auto pUnitSphere = Model::createFromFile(ModelUnitSphere, Model::LoadFlags::None);
	auto pUnitCylinder = Model::createFromFile(ModelUnitCylinder, Model::LoadFlags::None);
	auto pUnitCone = Model::createFromFile(ModelUnitCone, Model::LoadFlags::None);

	mpDbgScene = Scene::create();
	mpDbgSceneRenderer = SceneRenderer::create(mpDbgScene);

	for (int32_t y = 0; y < mHeight; ++y)
	{
		for (int32_t x = 0; x < mWidth; ++x)
		{
			auto pParticle = getParticle(x, y);

			pParticle->pDbgUnitSphere = Scene::ModelInstance::create(pUnitSphere, pParticle->mPosition, vec3(0.0f), vec3(0.001f));
			pParticle->pDbgUnitSphere->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
			mpDbgScene->addModelInstance(pParticle->pDbgUnitSphere);

			pParticle->pDbgUnitCylinder = Scene::ModelInstance::create(pUnitCylinder, pParticle->mPosition, vec3(0.0f), vec3(0.1f));
			pParticle->pDbgUnitCylinder->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
			mpDbgScene->addModelInstance(pParticle->pDbgUnitCylinder);

			// TODO: add cone to the end of cylinder
		}
	}
}

void init(const Model*)
{
	assert(!"Not implemented!");
}

ClothParticle* ClothPatch::getParticle(int32_t x, int32_t y)
{
	return &mpParticles[y * mWidth + x];
}

static constexpr struct
{
	float32 mTemperature, mDensity;
} AirPropertiesTable[] = {
	{ -25.0f, 1.4224f },
	{ -20.0f, 1.3943f },
	{ -15.0f, 1.3673f },
	{ -10.0f, 1.3413f },
	{ - 5.0f, 1.3163f },
	{   0.0f, 1.2922f },
	{   5.0f, 1.2690f },
	{  10.0f, 1.2466f },
	{  15.0f, 1.2250f },
	{  20.0f, 1.2041f },
	{  25.0f, 1.1839f },
	{  30.0f, 1.1644f },
	{  35.0f, 1.1455f },
};

static float32 GetAirDensity(float32 temperature)
{
	temperature = glm::clamp(
		temperature,
		AirPropertiesTable[0].mTemperature,
		AirPropertiesTable[arraysize(AirPropertiesTable) - 1].mTemperature
	);

	for (int32_t i = 1; i < arraysize(AirPropertiesTable); ++i)
	{
		const auto &temperature0 = AirPropertiesTable[i - 1].mTemperature;
		const auto &temperature1 = AirPropertiesTable[i].mTemperature;

		if (temperature <= temperature1)
		{
			return lerp(
				AirPropertiesTable[i - 1].mDensity,
				AirPropertiesTable[i].mDensity,
				(temperature - temperature0) / (temperature1 - temperature0)
			);
		}
	}

	return AirPropertiesTable[0].mDensity;
}

void ClothPatch::simulate(float deltaTime)
{
	const auto dragCoefficient = 0.47f;
	const auto area = 1.0f;
	const auto airDensity = GetAirDensity(mAirTemperature);

	const auto windAcceleration = float3(3.0f, 0.0f, 3.0f) / mMass;
	const auto gravityAcceleration = float3(0.0f, -9.8067f, 0.0f);

	for (int32_t y = 0; y < mHeight; ++y)
	{
		for (int32_t x = 0; x < mWidth; ++x)
		{
			auto particle = getParticle(x, y);
			if (particle->mbStationary) continue;

			const auto currPosition = particle->mPosition;
			const auto velocity = particle->mPosition - particle->mPrevPosition;

			const auto dragForce = 0.5f * airDensity * dragCoefficient * area * (velocity * velocity);
			const auto dragAcceleration = dragForce * mMass;

			const auto acceleration = gravityAcceleration + windAcceleration /*+ dragAcceleration*/;

			particle->mPosition += velocity + acceleration * deltaTime * deltaTime;
			particle->mPrevPosition = currPosition;
		}
	}

	const int32_t Offsets[][2] = {
		{  0,  1 },
		{  0, -1 },
		{  1,  0 },
		{ -1,  0 },

		{  1,  1 },
		{  1, -1 },
		{ -1,  1 },
		{ -1, -1 },

		{  0,  2 },
		{  0, -2 },
		{  2,  0 },
		{ -2,  0 },
	};

	for (int32_t y = 0; y < mHeight; ++y)
	{
		for (int32_t x = 0; x < mWidth; ++x)
		{
			auto first = getParticle(x, y);
			float3 adjust(0.0f);

			for (const auto &offset : Offsets)
			{
				if (x + offset[0] < 0 || x + offset[0] >= mWidth) continue;
				if (y + offset[1] < 0 || y + offset[1] >= mHeight) continue;

				auto second = getParticle(x + offset[0], y + offset[1]);
				auto distance = first->mPosition - second->mPosition;
				auto current = length(distance);
				auto target = sqrtf(float32_t(offset[0] * offset[0] + offset[1] * offset[1]));

				float weightFirst, weightSecond;
				if (first->mbStationary && second->mbStationary)
					weightFirst = 0.0f, weightSecond = 0.0f;
				else if (!first->mbStationary && second->mbStationary)
					weightFirst = 1.0f, weightSecond = 0.0f;
				else if (first->mbStationary && !second->mbStationary)
					weightFirst = 0.0f, weightSecond = 1.0f;
				else
					weightFirst = 0.5f, weightSecond = 0.5f;

				distance *= ((target - current) / current);
				first->mPosition += distance * weightFirst;
				second->mPosition -= distance * weightSecond;
			}
		}
	}
}

void ClothPatch::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	const auto currTime = pSample->getCurrentTime();
	auto *pRenderContext = pSample->getRenderContext();
	const auto &pTargetFbo = pSample->getCurrentFbo();

	float32_t *pPositions = reinterpret_cast<float32_t*>(mpVBPositions->map(Buffer::MapType::WriteDiscard));
	float32_t *pNormals = reinterpret_cast<float32_t*>(mpVBNormals->map(Buffer::MapType::WriteDiscard));
	for (int32_t y = 0; y < mHeight; ++y)
	{
		for (int32_t x = 0; x < mWidth; ++x)
		{
			const auto *particle = getParticle(x, y);
			auto *position = pPositions + (y * mWidth + x) * 3;

			position[0] = particle->mPosition[0];
			position[1] = particle->mPosition[1];
			position[2] = particle->mPosition[2];

			// TODO: Fit a plane to the particles' position and use its normal

			// Normal sampling pattern:
			// *  *  *
			//   /|\
            //  /4|3\
            // *__*__*
			//  \1|2/
			//   \|/
			// *  *  *

			const int32_t NormalSamples[][3][2] =
			{
				{ { -1,  0 }, {  0,  0 }, {  0,  1 } }, // 1
				{ {  0,  1 }, {  0,  0 }, {  1,  0 } }, // 2
				{ {  1,  0 }, {  0,  0 }, {  0, -1 } }, // 3
				{ {  0, -1 }, {  0,  0 }, { -1,  0 } }, // 4
			};

			const auto &GetParticleSample = [&](const int32_t(&offset)[2]) -> float3 {

				int32_t x1 = glm::clamp(x + offset[0], 0, mWidth);
				int32_t y1 = glm::clamp(y + offset[1], 0, mHeight);
				return getParticle(x1, y1)->mPosition;
			};
			const auto &GetPlaneNormal = [](const float3 &A, const float3 &B, const float3 &C) -> float3 {
				const auto result = normalize(cross(C - B, A - B));
				return any(isnan(result)) ? float3(0.0f) : result;
			};

			float3 accumulator(0);
			for (const auto &sample : NormalSamples)
			{
				accumulator += GetPlaneNormal(
					GetParticleSample(sample[0]), GetParticleSample(sample[1]), GetParticleSample(sample[2])
				);
			}

			accumulator = normalize(accumulator);
			auto outNormal = pNormals + (y * mWidth + x) * 3;
			outNormal[0] = accumulator[0];
			outNormal[1] = accumulator[1];
			outNormal[2] = accumulator[2];

			if (mbShowParticles)
			{
				const auto velocity = particle->mPosition - particle->mPrevPosition;
				particle->pDbgUnitSphere->setTranslation(particle->mPosition, true);
				particle->pDbgUnitCylinder->setTranslation(particle->mPosition, true);
				particle->pDbgUnitCylinder->setRotation(getYawPitchRoll(velocity));
				particle->pDbgUnitCylinder->setScaling(vec3(0.1, 0.1, length(velocity)));
			}
		}
	}
	mpVBPositions->unmap();
	mpVBNormals->unmap();

	const auto &GetTriangles = [&](int32_t x, int32_t y, vec3 triangle0[3], vec3 triangle1[3]) -> void
	{
		const auto &GetVertex = [&](int32_t x, int32_t y) -> vec3
		{
			return getParticle(x, y)->mPosition;
		};

		triangle0[0] = GetVertex(x, y + 1);
		triangle0[1] = GetVertex(x + 1, y);
		triangle0[2] = GetVertex(x, y);

		triangle1[0] = GetVertex(x + 1, y + 1);
		triangle1[1] = GetVertex(x + 1, y);
		triangle1[2] = GetVertex(x, y + 1);
	};

	auto *pColors = reinterpret_cast<vec3*>(mpVBColors->map(Buffer::MapType::WriteDiscard));
	std::fill_n(pColors, mHeight * mWidth, vec3(0.0f, 1.0f, 0.0f));

	const auto &HandleCollision = [&](int32_t x, int32_t y) -> void
	{
		pColors[y * mWidth + x] = vec3(1.0f, 0.0f, 0.0f);
		// getParticle(x, y)->mPosition = getParticle(x, y)->mPrevPosition;
	};

	for (const auto &triangle : triangles)
	{
		const vec3 vertices[3] = {
			mpParticles[triangle.v[0]].mPosition,
			mpParticles[triangle.v[1]].mPosition,
			mpParticles[triangle.v[2]].mPosition
		};
	}

	for (int32_t y0 = 0; y0 < mHeight - 1; ++y0)
	{
		for (int32_t x0 = 0; x0 < mWidth - 1; ++x0)
		{
			vec3 triangle0_0[3], triangle0_1[3];
			GetTriangles(x0, y0, triangle0_0, triangle0_1);

			for (int32_t y1 = y0; y1 < mHeight - 1; ++y1)
			{
				for (int32_t x1 = y1 == y0 ? x0 + 1 : 0; x1 < mWidth - 1; ++x1)
				{
					if (std::abs(y0 - y1) <= 1 && std::abs(x0 - x1) <= 1)
						continue;

					vec3 triangle1_0[3], triangle1_1[3];
					GetTriangles(x1, y1, triangle1_0, triangle1_1);

					vec3 responses[4];

					bool collisions[4] =
					{
						intersectTriangles(triangle0_0, triangle1_0, responses[0]),
						intersectTriangles(triangle0_0, triangle1_1, responses[1]),
						intersectTriangles(triangle0_1, triangle1_0, responses[2]),
						intersectTriangles(triangle0_1, triangle1_1, responses[3])
					};

					if (collisions[0] || collisions[1])
					{
						HandleCollision(x0, y0 + 1);
						HandleCollision(x0 + 1, y0);
						HandleCollision(x0, y0);
					}
					if (collisions[2] || collisions[3])
					{
						HandleCollision(x0 + 1, y0 + 1);
						HandleCollision(x0 + 1, y0);
						HandleCollision(x0, y0 + 1);
					}
					if (collisions[0] || collisions[2])
					{
						HandleCollision(x1, y1 + 1);
						HandleCollision(x1 + 1, y1);
						HandleCollision(x1, y1);
					}
					if (collisions[1] || collisions[3])
					{
						HandleCollision(x1 + 1, y1 + 1);
						HandleCollision(x1 + 1, y1);
						HandleCollision(x1, y1 + 1);
					}
				}
			}
		}
	}
	mpVBColors->unmap();

	if (mbShowWireframe)
		mpClothState->setRasterizerState(mpRasterizeWireframe);
	else
		mpClothState->setRasterizerState(mpRasterizeNormal);

	pClothSample->mpCamera->setIntoConstantBuffer(mpClothVars->getConstantBuffer("InternalPerFrameCB").get(), "gCamera");
	pClothSample->mpDirLight->setIntoProgramVars(mpClothVars.get(), mpClothVars["PerFrameCB"].get(), "gDirLight");

	pRenderContext->setGraphicsVars(mpClothVars);
	pRenderContext->setGraphicsState(mpClothState);

	mpClothState->pushFbo(pTargetFbo);
	pRenderContext->drawIndexed(triangles.size() * 3, 0, 0);
	mpClothState->popFbo();

	if (mbShowParticles)
	{
		mpDbgSceneRenderer->update(currTime);
		pClothSample->mpDirLight->setIntoProgramVars(mpModelVars.get(), mpModelVars["PerFrameCB"].get(), "gDirLight");

		pRenderContext->setGraphicsVars(mpModelVars);
		pRenderContext->setGraphicsState(mpModelState);

		mpModelState->pushFbo(pTargetFbo);
		mpDbgSceneRenderer->renderScene(pRenderContext, pClothSample->mpCamera.get());
		mpModelState->popFbo();
	}
}

void ClothPatch::onGuiRender(ClothSample*, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

	if (pGui->beginGroup(mName))
	{
		pGui->addCheckBox("Show Particles", mbShowParticles);
		pGui->addCheckBox("Show Wireframe", mbShowWireframe);
		pGui->endGroup();
	}
}

bool ClothPatch::onMouseEvent(ClothSample *pClothSample, SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
	const auto &cameraData = pClothSample->mpCamera->getData();
	const auto &GetMouseRay = [&](vec3 &rayOrigin, vec3 &rayDirection) -> void
	{
		rayOrigin = cameraData.posW;
		rayDirection = cameraData.cameraW;
		rayDirection += cameraData.cameraU * (mouseEvent.pos.x - 0.5f) * 2.0f;
		rayDirection += cameraData.cameraV * (0.5f - mouseEvent.pos.y) * 2.0f;
		rayDirection = normalize(rayDirection);
	};

	switch (mouseEvent.type)
	{
	case MouseEvent::Type::LeftButtonDown:
	{
		if (mHeldParticle == nullptr)
		{
			vec3 rayOrigin, rayDirection;
			GetMouseRay(rayOrigin, rayDirection);

			auto closestHit = std::numeric_limits<float32>::max();
			
			for (const auto &triangle : triangles)
			{
				const vec3 vertices[3] = {
					mpParticles[triangle.v[0]].mPosition,
					mpParticles[triangle.v[1]].mPosition,
					mpParticles[triangle.v[2]].mPosition
				};

				vec3 intersection;
				if (intersectTriangle(vertices, rayOrigin, rayDirection, intersection))
				{
					auto currentHit = 0.0f;
					for (int32 i = 0; i < 3; ++i)
					{
						if ((currentHit = length(intersection - vertices[i])) < closestHit)
						{
							mHeldParticle = &mpParticles[triangle.v[i]];
							closestHit = currentHit;
						}
					}
				}
			}

			if (mHeldParticle != nullptr)
			{
				mHeldParticle->mbStationary = true;
				return false;
			}
		}
	} break;
	case MouseEvent::Type::LeftButtonUp:
	{
		if (mHeldParticle != nullptr)
		{
			mHeldParticle->mbStationary = false;
			mHeldParticle = nullptr;
			return false;
		}

	} break;
	}

	if (mHeldParticle != nullptr)
	{
		vec3 rayOrigin, rayDirection;
		GetMouseRay(rayOrigin, rayDirection);

		const auto UxV = cross(cameraData.cameraU, cameraData.cameraV);
		const auto t = dot(rayOrigin - mHeldParticle->mPosition, UxV) / dot(-rayDirection, UxV);
		mHeldParticle->mPosition = rayOrigin + t * rayDirection;
	}
}