#include "External/HalfEdge/trimesh.h"

#include "TriangleIntersection.h"
#include "ClothSample.h"

using float1 = float32_t;
using float3 = f32vec3;

struct ClothParticle
{
	float3 mPosition;
	float3 mPrevPosition;
	float3 mAcceleration;
	bool mbStationary;
};

struct ClothPatch : public ClothModel
{
	ClothParticle *mpParticles;
	int32_t mWidth;
	int32_t mHeight;
	float32 mMass;

	ClothParticle *mHeldParticle;
	float32 mAirTemperature;

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
	vec3 getVertexPosition(int32 vertexIndex) const final override;

	EType getType() const final override
	{
		return EType::ParticleSpringModel;
	}

	void init(int32_t width, int32_t height) final override;
	void init(const Model *model) final override;
	void simulate(float deltaTime) final override;
	void render(ClothSample*, SampleCallbacks*) final override;

	void onGuiRender(ClothSample*, SampleCallbacks*) final override;
	bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) final override;
};

void ClothPatch::init(int32_t width, int32_t height)
{
	std::vector<float3> positions;
	const auto vertexCount = createRectMesh(width, height, positions, mTriangles, mMesh);

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

	getParticle(0, 0)->mbStationary = true;
	getParticle(width - 1, 0)->mbStationary = true;
	
	sharedInit();
}

void ClothPatch::init(const Model*)
{
	assert(!"Not implemented!");
}

ClothParticle* ClothPatch::getParticle(int32_t x, int32_t y)
{
	return &mpParticles[y * mWidth + x];
}

vec3 ClothPatch::getVertexPosition(int32 vertexIndex) const
{
	return mpParticles[vertexIndex].mPosition;
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

	const auto windAcceleration = float3(3.0f, 0.0f, 3.0f) / mMass; // TODO: Take into account mMaterialDensity
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
			const auto dragAcceleration = dragForce * mMass; // TODO: Take into account mMaterialDensity

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
	ClothModel::render(pClothSample, pSample);

	const auto currTime = pSample->getCurrentTime();
	auto *pRenderContext = pSample->getRenderContext();
	const auto &pTargetFbo = pSample->getCurrentFbo();

	if (mbShowParticles)
	{
		auto dbgScene = Scene::create();
		auto dbgSceneRenderer = SceneRenderer::create(dbgScene);

		for (int32 i = 0; i < mWidth * mHeight; ++i)
		{
			const auto &particle = mpParticles[i];
			const auto velocity = particle.mPosition - particle.mPrevPosition;

			auto unitCylinder = Scene::ModelInstance::create(
				pClothSample->mpDbgUnitCylinder,
				particle.mPosition,
				getYawPitchRoll(velocity),
				vec3(0.1, 0.1, length(velocity))
			);
			unitCylinder->move(vec3(0.0f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 1.0f, 0.0f));
			dbgScene->addModelInstance(unitCylinder);
		}

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

void ClothPatch::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

	if (pGui->beginGroup(mName))
	{
		ClothModel::onGuiRender(pClothSample, pSample);
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
			
			for (const auto &triangle : mTriangles)
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
				return true;
			}
		}
	} break;
	case MouseEvent::Type::LeftButtonUp:
	{
		if (mHeldParticle != nullptr)
		{
			mHeldParticle->mbStationary = false;
			mHeldParticle = nullptr;
			return true;
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
		
		return true;
	}

	return false;
}

template <>
ClothModel* ClothModel::createClothModel<ClothModel::ParticleSpringModel>()
{
	return new ::ClothPatch();
}