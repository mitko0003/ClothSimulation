#include "External/HalfEdge/trimesh.h"

#include "TriangleIntersection.h"
#include "ClothSample.h"
#include "glm/gtc/random.hpp"
#pragma optimize("", off)
using float1 = float32_t;
using float3 = f32vec3;

struct IBaseParticle
{
    float3 mPosition;
    float3 mSurfaceNormal;
    float mSurfaceArea;
    bool mbStationary;

    IBaseParticle() :
        mPosition(0.0f), mbStationary(false) {}
    virtual ~IBaseParticle() {}

    virtual void simulate(float deltaTime, float3 acceleration) = 0;
    virtual float3 getVelocity(float deltaTime) const = 0;

    virtual void setPosition(const float3 &position)
    {
        mPosition = position;
    }
};

struct TEulerParticle : public IBaseParticle
{
    float3 mVelocity;

    TEulerParticle() :
        mVelocity(0.0f) {}

    void simulate(float deltaTime, float3 acceleration) final override
    {
        mVelocity += acceleration * deltaTime;
        mPosition += mVelocity * deltaTime;
    }

    float3 getVelocity(float) const final override
    {
        return mVelocity;
    }
};

struct IVerletWithoutVelocityParticle : public IBaseParticle
{
    float3 mPrevPosition;

    IVerletWithoutVelocityParticle() :
        mPrevPosition(mPosition) {}

    float3 getVelocity(float deltaTime) const final override
    {
        return (mPosition - mPrevPosition) / deltaTime;
    }

    virtual void setPosition(const float3 &position) final override
    {
        mPrevPosition = mPosition = position;
    }
};

struct TStormerVerletParticle : public IVerletWithoutVelocityParticle
{
    using IVerletWithoutVelocityParticle::IVerletWithoutVelocityParticle;

    void simulate(float deltaTime, float3 acceleration) final override
    {
        const auto prevPosition = mPosition;
        mPosition += mPosition - mPrevPosition + acceleration * Square(deltaTime);
        mPrevPosition = prevPosition;
    }
};

struct TVerletBasicNonConstTimeDeltaParticle : public IVerletWithoutVelocityParticle
{
    using IVerletWithoutVelocityParticle::IVerletWithoutVelocityParticle;

    float32 mPrevDeltaTime = 0.01f;

    void simulate(float deltaTime, float3 acceleration) final override
    {
        const auto prevPosition = mPosition;
        mPosition = mPosition + (mPosition - mPrevPosition) * deltaTime / mPrevDeltaTime + acceleration * Square(deltaTime);
        mPrevPosition = prevPosition;
        mPrevDeltaTime = deltaTime;
    }
};

struct TVerletImprovedNonConstTimeDeltaParticle : public IVerletWithoutVelocityParticle
{
    using IVerletWithoutVelocityParticle::IVerletWithoutVelocityParticle;

    float32 mPrevDeltaTime = 0.01f;

    void simulate(float deltaTime, float3 acceleration) final override
    {
        const auto prevPosition = mPosition;
        mPosition = mPosition + (mPosition - mPrevPosition) * (deltaTime / mPrevDeltaTime) + acceleration * ((deltaTime + mPrevDeltaTime) * deltaTime * 0.5f);
        mPrevPosition = prevPosition;
        mPrevDeltaTime = deltaTime;
    }
};

struct TVelocityVerletParticle : public IBaseParticle
{
    float3 mVelocity;
    float3 mAcceleration;

    TVelocityVerletParticle() :
        mVelocity(0.0f), mAcceleration(0.0f) {}

    void simulate(float deltaTime, float3 acceleration) final override
    {
        const auto halfStepVelocity = getVelocity(deltaTime);
        mPosition = mPosition + halfStepVelocity * deltaTime;
        mVelocity = halfStepVelocity + acceleration * (deltaTime * 0.5f);
        mAcceleration = acceleration;
    }

    float3 getVelocity(float deltaTime) const final override
    {
        return mVelocity + mAcceleration * (deltaTime * 0.5f);
    }
};

enum ENumericalMethod
{
    ExplicitEuler,
    RungeKutta,
	StormerVerlet,
    VelocityVerlet,
	VerletBasicNonConstTimeDelta,
	VerletImprovedNonConstTimeDelta,

	NumericalMethodCount
};

struct
{
	const char *mpName;
    int32 mParticleSize;
	bool mbConstantTimeDelta;
} NumericalMethodProperties[NumericalMethodCount] = {
	{ "Explicit Euler Method", sizeof(TEulerParticle), true },
	{ "Runge-Kutta Method", sizeof(TEulerParticle), true },
	{ "Stormer-Verlet Method", sizeof(TStormerVerletParticle), true },
	{ "Velocity Verlet Method", sizeof(TVelocityVerletParticle), true },
	{ "Verlet non-const time delta", sizeof(TVerletBasicNonConstTimeDeltaParticle), false },
	{ "Verlet non-const time delta improved", sizeof(TVerletImprovedNonConstTimeDeltaParticle), false },
};

enum ESpringType
{
	SpringTypeStructural,
	SpringTypeShear,
	SpringTypeFlexion,

	SpringTypeCount
};

static constexpr struct
{
	int32 mOffset[2];
	ESpringType mType;
} Springs[] =
{
	{ {  0,  1 }, SpringTypeStructural },
	{ {  0, -1 }, SpringTypeStructural },
	{ {  1,  0 }, SpringTypeStructural },
	{ { -1,  0 }, SpringTypeStructural },

	{ {  1,  1 }, SpringTypeShear },
	{ {  1, -1 }, SpringTypeShear },
	{ { -1,  1 }, SpringTypeShear },
	{ { -1, -1 }, SpringTypeShear },

	{ {  0,  2 }, SpringTypeFlexion },
	{ {  0, -2 }, SpringTypeFlexion },
	{ {  2,  0 }, SpringTypeFlexion },
	{ { -2,  0 }, SpringTypeFlexion },
};

struct ParticleSpringModel final : public ClothModel
{
	uint8 *mpParticlesMemory = nullptr;
	vec2 mSize = vec2(0, 0);
	ivec2 mTessellation = ivec2(0, 0);
    float mTimeRemainder = 0.0f;

    IBaseParticle *mHeldParticle = nullptr;

	struct
	{
		float32 airTemperature = 20.0f;
		float32 timeStep = 0.1f;
		float32 particleMass = 0.1f;
        float32 springStiffness = 0.0f;
        float32 dampingCoefficient = 0.1f;
		ENumericalMethod numericalMethod = StormerVerlet;
		bool useSprings[SpringTypeCount];
	} mUserParams;

	ParticleSpringModel()
	{
		mName = "Particle Spring Model";
		std::fill_n(mUserParams.useSprings, arraysize(mUserParams.useSprings), true);
	}

	~ParticleSpringModel()
	{
		delete[] mpParticlesMemory;
	}

    IBaseParticle* getParticle(int32 index) {
        const auto particleSize = NumericalMethodProperties[mUserParams.numericalMethod].mParticleSize;
        return reinterpret_cast<IBaseParticle*>(mpParticlesMemory + index * particleSize);
    }
    const IBaseParticle* getParticle(int32 index) const {
        return const_cast<ParticleSpringModel*>(this)->getParticle(index);
    }
    IBaseParticle* getParticle(int32 x, int32 y) {
        return getParticle(y * mTessellation.x + x);
    }
    const IBaseParticle* getParticle(int32 x, int32 y) const {
        return const_cast<ParticleSpringModel*>(this)->getParticle(x, y);
    }

    void updateSurfaceNormals();
    float32_t ParticleSpringModel::computeVoronoiArea(int32 particle) const;

    vec3 getVertexNormal(int32 index) const override {
        return getParticle(index)->mSurfaceNormal;
    }
	EType getType() const final override {
		return EType::ParticleSpringModel;
	}

	void init(const vec2 &size, const ivec2 &tessellation) override;
	void init(const Model *model) override;
	void simulate(ClothSample*, float deltaTime) override;
	void render(ClothSample*, SampleCallbacks*) override;
    vec3 getVertexPosition(int32 vertexIndex) const override;

	void onGuiRender(ClothSample*, SampleCallbacks*) override;

    void testSelection(SelectionQuery&) override;
    void makeSelection(ClothSample*, SelectionQuery&) override;
    void loseSelection(ClothSample*) override;

	bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) override;
};

static uint8* createParticles(ENumericalMethod numericalMethod, int32 count, std::function<void(int32, IBaseParticle*)> initParticle)
{
    const auto particleSize = NumericalMethodProperties[numericalMethod].mParticleSize;
    auto *particlesMemory = new uint8[count * particleSize];

    for (int32_t i = 0; i < count; ++i)
    {
        auto *newParticle = reinterpret_cast<IBaseParticle*>(particlesMemory + i * particleSize);
        switch (numericalMethod)
        {
        case ExplicitEuler:
        case RungeKutta: new (newParticle) TEulerParticle(); break;
        case VelocityVerlet: new (newParticle) TVelocityVerletParticle(); break;
        case StormerVerlet: new (newParticle) TStormerVerletParticle(); break;
        case VerletBasicNonConstTimeDelta: new (newParticle) TVerletBasicNonConstTimeDeltaParticle(); break;
        case VerletImprovedNonConstTimeDelta: new (newParticle) TVerletImprovedNonConstTimeDeltaParticle(); break;
        default:
            should_not_get_here();
        }
        initParticle(i, newParticle);
    }

    return particlesMemory;
}

void ParticleSpringModel::updateSurfaceNormals()
{
    const auto vertexCount = int32(mMesh.get_num_vertices());
    for (int32 vertex = 0; vertex < vertexCount; ++vertex)
        getParticle(vertex)->mSurfaceNormal = computeVertexNormal(vertex);
}

float32_t ParticleSpringModel::computeVoronoiArea(int32 particle) const
{
    std::vector<float3> vertices;
    auto faces = mMesh.vertex_face_neighbors(particle);

    auto area = 0.0f;
    for (const auto &face : faces)
    {
        const auto &triangle = mTriangles[face];

        int32 rotate = -1;
        while (triangle[++rotate] != particle);

        const auto xi = getParticle(triangle[(rotate + 0) % 3])->mPosition;
        const auto xj = getParticle(triangle[(rotate + 1) % 3])->mPosition;
        const auto xk = getParticle(triangle[(rotate + 2) % 3])->mPosition;

        const auto xij = xj - xi;
        const auto xik = xk - xi;
        const auto xjk = xk - xj;

        const auto dir_ij = normalize(xij);
        const auto dir_ik = normalize(xik);
        const auto dir_jk = normalize(xjk);

        const auto cosPhi_i = dot(dir_ij, dir_ik);
        const auto cosPhi_j = -dot(dir_jk, dir_ij);
        const auto cosPhi_k = dot(dir_jk, dir_ik);

        const auto isTriangleObtuse = cosPhi_i < 0.0f || cosPhi_j < 0.0f || cosPhi_k < 0.0f;
        if (!isTriangleObtuse)
        {
            const auto cotgPhi_j = cosPhi_j / sqrtf(1.0f - cosPhi_j * cosPhi_j);
            const auto cotgPhi_k = cosPhi_k / sqrtf(1.0f - cosPhi_k * cosPhi_k);

            area += (dot(xik, xik) * cotgPhi_j + dot(xij, xij) * cotgPhi_k) / 8.0f;
        }
        else
        {
            const auto areaTriangle = length(cross(xij, xik)) / 2.0f;

            const auto isXiObtuse = cosPhi_i < 0.0f;
            if (isXiObtuse)
                area += areaTriangle / 2.0f;
            else
                area += areaTriangle / 4.0f;
        }
    }
    return area;
}

void ParticleSpringModel::init(const vec2 &size, const ivec2 &tessellation)
{
	std::vector<float3> positions;
	const auto vertexCount = createRectMesh(size, tessellation, positions, mTriangles, mMesh);

	mSize = size;
	mTessellation = tessellation;
    mpParticlesMemory = createParticles(mUserParams.numericalMethod, vertexCount,
        [&](int32 index, IBaseParticle *particle) {
            particle->setPosition(positions[index]);
        }
    );

	getParticle(0, 0)->mbStationary = true;
	getParticle(tessellation.x - 1, 0)->mbStationary = true;
    updateSurfaceNormals();

    for (int32 vertex = 0; vertex < vertexCount; ++vertex)
        getParticle(vertex)->mSurfaceArea = computeVoronoiArea(vertex);

	sharedInit();
}

void ParticleSpringModel::init(const Model*)
{
	assert(!"Not implemented!");
}

vec3 ParticleSpringModel::getVertexPosition(int32 vertexIndex) const
{
	return getParticle(vertexIndex)->mPosition;
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

void ParticleSpringModel::simulate(ClothSample *pClothSample, float deltaTime)
{
    mTimeRemainder += deltaTime;
	std::function<bool()> TimeIteration;
	if (NumericalMethodProperties[mUserParams.numericalMethod].mbConstantTimeDelta)
	{
		TimeIteration = [&]() -> bool
		{
			if (mTimeRemainder >= mUserParams.timeStep)
			{
                deltaTime = mUserParams.timeStep;
                mTimeRemainder -= mUserParams.timeStep;
				return true;
			}
			return false;
		};
	}
	else
	{
		TimeIteration = [&]() -> bool
		{
			if (mTimeRemainder > 0.0001f)
			{
                deltaTime = min(mTimeRemainder, mUserParams.timeStep);
                mTimeRemainder -= deltaTime;
				return true;
			}
			return false;
		};
	}

	const auto dragCoefficient = 0.47f;
	const auto area = mSize.x * mSize.y / (mTessellation.x * mTessellation.y);
	const auto ro = GetAirDensity(mUserParams.airTemperature);
    const auto gridStep = mSize / vec2(mTessellation - 1);
	const auto gravityAcceleration = float3(0.0f, -9.8067f, 0.0f);
    const auto windVelocity = pClothSample->getWindVelocity();
    const auto Cd = ClothModel::mUserParams.windDragCoefficient;
    const auto Cl = ClothModel::mUserParams.windLiftCoefficient;

    while (TimeIteration())
    {
		for (int32_t y = 0; y < mTessellation.y; ++y)
		{
			for (int32_t x = 0; x < mTessellation.x; ++x)
			{
				auto particle = getParticle(x, y);
				if (particle->mbStationary) continue;

                vec3 springForce(0.0f);
                if (mUserParams.springStiffness > 0.0f)
                {
                    for (const auto &spring : Springs)
                    {
                        if (spring.mType == SpringTypeFlexion)
                            continue;

                        const auto &offset = spring.mOffset;
                        if (x + offset[0] < 0 || x + offset[0] >= mTessellation.x) continue;
                        if (y + offset[1] < 0 || y + offset[1] >= mTessellation.y) continue;

                        auto second = getParticle(x + offset[0], y + offset[1]);
                        auto distance = particle->mPosition - second->mPosition;
                        auto current = length(distance);
                        auto target = sqrtf(Square(offset[0] * gridStep.x) + Square(offset[1] * gridStep.y));

                        springForce -= distance * (mUserParams.springStiffness * (1.0f - target / current));
                    }
                }

                const auto v_p = particle->getVelocity(deltaTime);
                const auto v_r = windVelocity - v_p;
                const auto &A = particle->mSurfaceArea;
                auto n = particle->mSurfaceNormal;
                if (dot(n, v_r) > 0.0f) n *= -1; // Flip the surface normal based on which side the wind blows on.

                // Equation from: "Simulating Wind Effects on Cloth and Hair in Disney's Frozen"
                const auto windForce = 0.5f * ro * A * ((Cd - Cl) * dot(v_r, n) * v_r + Cl * length(v_p) * n);
				const auto windAcceleration = windForce / mUserParams.particleMass;
                
                const auto dampingForce = -mUserParams.dampingCoefficient * v_p;

				const auto acceleration = gravityAcceleration + windAcceleration + dampingForce;
                particle->simulate(deltaTime, acceleration);

				if (ClothModel::mUserParams.bShowWindEffect)
					pClothSample->drawVector(windAcceleration, particle->mPosition);
			}
		}

        for (int32_t y = 0; y < mTessellation.y; ++y)
        {
            for (int32_t x = 0; x < mTessellation.x; ++x)
            {
                auto first = getParticle(x, y);
                float3 adjust(0.0f);

                for (const auto &spring : Springs)
                {
                    if (!mUserParams.useSprings[spring.mType])
                        continue;

                    const auto &offset = spring.mOffset;
                    if (x + offset[0] < 0 || x + offset[0] >= mTessellation.x) continue;
                    if (y + offset[1] < 0 || y + offset[1] >= mTessellation.y) continue;

                    auto second = getParticle(x + offset[0], y + offset[1]);
                    auto distance = second->mPosition - first->mPosition;
                    auto current = length(distance);
                    auto target = sqrtf(Square(offset[0] * gridStep.x) + Square(offset[1] * gridStep.y));

                    if (current > target * 1.1f)
                        target = target * 1.1f;
                    else if (current < target * 0.9f)
                        target = target * 0.9f;
                    else
                        continue;

                    float weightFirst, weightSecond;
                    if (first->mbStationary && second->mbStationary)
                        weightFirst = 0.0f, weightSecond = 0.0f;
                    else if (!first->mbStationary && second->mbStationary)
                        weightFirst = 1.0f, weightSecond = 0.0f;
                    else if (first->mbStationary && !second->mbStationary)
                        weightFirst = 0.0f, weightSecond = 1.0f;
                    else
                        weightFirst = 0.5f, weightSecond = 0.5f;

                    distance *= ((current - target) / current);
                    first->mPosition += distance * weightFirst;
                    second->mPosition -= distance * weightSecond;
                }
            }
        }

        for (const auto &triangle : mTriangles)
        {
            vec3 vertices[3];
            getTriangleVertices(triangle, vertices);

            for (const auto &object : pClothSample->mPhysicsObjects)
            {
                const auto response = object->intersectionResponse(vertices);
                if (response != vec3(0.0f))
                {
                    getParticle(triangle.i)->mPosition += response;
                    getParticle(triangle.j)->mPosition += response;
                    getParticle(triangle.k)->mPosition += response;
                }
            }
        }

        updateSurfaceNormals();
    }
}

void ParticleSpringModel::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	ClothModel::render(pClothSample, pSample);

	const auto currTime = pSample->getCurrentTime();
	auto *pRenderContext = pSample->getRenderContext();
	const auto &pTargetFbo = pSample->getCurrentFbo();
}

void ParticleSpringModel::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

	if (pGui->beginGroup(mName, true))
	{

		Gui::DropdownList numericalMethodDropdown;
		for (uint32 i = 0; i < NumericalMethodCount; ++i)
			numericalMethodDropdown.push_back({ i, NumericalMethodProperties[i].mpName });

        uint32 newNumericalMethod = mUserParams.numericalMethod;
        if (pGui->addDropdown("Numerical Method", numericalMethodDropdown, newNumericalMethod))
        {
            auto newParticles = createParticles(ENumericalMethod(newNumericalMethod), mTessellation.x * mTessellation.y,
                [&](int32 index, IBaseParticle *particle) {
                    const auto *oldParticle = getParticle(index);
                    particle->setPosition(oldParticle->mPosition);
                    particle->mbStationary = oldParticle->mbStationary;
                }
            );
            delete[] mpParticlesMemory;
            mpParticlesMemory = newParticles;
            mUserParams.numericalMethod = ENumericalMethod(newNumericalMethod);
        }

		if (NumericalMethodProperties[mUserParams.numericalMethod].mbConstantTimeDelta)
            pClothSample->addFloatSlider("Const Time Step", mUserParams.timeStep, 0.0001f, 0.1f, false, "%.4f");
        else
            pClothSample->addFloatSlider("Max Time Step", mUserParams.timeStep, 0.0001f, 0.1f, false, "%.4f");
        pGui->addFloatSlider("Damping Coefficient", mUserParams.dampingCoefficient, 0.0f, 1.0f);
		pGui->addFloatSlider("Air Temperature", mUserParams.airTemperature, -25.0f, 35.0f);
		pGui->addFloatSlider("Particle Mass", mUserParams.particleMass, 0.0001f, 1.0f);
        pGui->addFloatSlider("Spring Stiffness", mUserParams.springStiffness, 0.0f, 10.0f);
		
		pGui->addCheckBox("Use Structural Springs", mUserParams.useSprings[SpringTypeStructural]);
		pGui->addCheckBox("Use Shear Springs", mUserParams.useSprings[SpringTypeShear]);
		pGui->addCheckBox("Use Flexion Springs", mUserParams.useSprings[SpringTypeFlexion]);

		ClothModel::onGuiRender(pClothSample, pSample);
		pGui->endGroup();
	}
}

struct SelectionQueryCache
{
    IBaseParticle *selectedParticle = nullptr;
};

void ParticleSpringModel::testSelection(SelectionQuery &query)
{
    auto closestHit = query.closestHit;
    SelectionQueryCache localSelectionCache;

    for (const auto &triangle : mTriangles)
    {
        const vec3 vertices[3] = {
            getParticle(triangle.v[0])->mPosition,
            getParticle(triangle.v[1])->mPosition,
            getParticle(triangle.v[2])->mPosition
        };

        vec3 intersection;
        if (intersectTriangle(vertices, query.ray.origin, query.ray.direction, intersection))
        {
            IBaseParticle *hitParticle = nullptr;
            auto minDistance = std::numeric_limits<float32>::max();

            for (int32 i = 0; i < 3; ++i)
            {
                auto currDistance = length(intersection - vertices[i]);
                if (currDistance < minDistance)
                {
                    hitParticle = getParticle(triangle.v[i]);
                    minDistance = currDistance;
                }
            }

            auto currentHit = length(query.ray.origin - hitParticle->mPosition);
            if (currentHit < closestHit)
            {
                localSelectionCache.selectedParticle = hitParticle;
                closestHit = currentHit;
            }
        }
    }

    if (localSelectionCache.selectedParticle != nullptr)
    {
        query.closestHit = closestHit;
        query.closestObject = static_cast<SelectableObject*>(this);

        auto *selectionCache = reinterpret_cast<SelectionQueryCache*>(query.cache);
        *selectionCache = localSelectionCache;
    }
}

void ParticleSpringModel::makeSelection(ClothSample*, SelectionQuery &query)
{
    auto *selectionCache = reinterpret_cast<SelectionQueryCache*>(query.cache);
    mHeldParticle = selectionCache->selectedParticle;
    mHeldParticle->mbStationary = true;
}

void ParticleSpringModel::loseSelection(ClothSample*)
{
    mHeldParticle->mbStationary = false;
    mHeldParticle = nullptr;
}

bool ParticleSpringModel::onMouseEvent(ClothSample *pClothSample, SampleCallbacks* pSample, const MouseEvent& mouseEvent)
{
	const auto &cameraData = pClothSample->mpCamera->getData();

	if (mHeldParticle != nullptr)
	{
		vec3 rayOrigin, rayDirection;
		GetMouseRay(mouseEvent, cameraData, rayOrigin, rayDirection);

		const auto UxV = cross(cameraData.cameraU, cameraData.cameraV);
		const auto t = dot(mHeldParticle->mPosition - rayOrigin, UxV) / dot(rayDirection, UxV);
		mHeldParticle->mPosition = rayOrigin + t * rayDirection;
		
		return true;
	}

	return false;
}

template <>
auto ClothModel::createClothModel<ClothModel::ParticleSpringModel>(const ClothModel *baseParent) -> SharedPtr
{
	auto result = new ::ParticleSpringModel();
	if (baseParent != nullptr)
	{
		const auto *parent = static_cast<const ::ParticleSpringModel*>(baseParent);
		result->mUserParams = parent->mUserParams;
	}

	return SharedPtr(result);
}
