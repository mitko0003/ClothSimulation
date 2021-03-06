#include "Eigen/Eigen"
#include "External/HalfEdge/trimesh.h"

#include "ClothSample.h"

static constexpr float64_t ShearModulus(float64_t E45, float64_t V45)
{
	return E45 * 0.5 / (1.0 + V45);
}

struct TMaterialProperties // For orthotropic materials
{
    const char *Name;
    float32_t Density;   // Density in kg/m^2
	float64_t Ex, Ey;    // Young’s modulus in KPa
	float64_t Vxy, Vyx;  // Poisson’s ratio
	float64_t Es;        // Shear modulus in KPa

	Eigen::Matrix<float32_t, 3, 3> getHookeLawMatrix() const
	{
		return Eigen::Matrix<float32_t, 3, 3>({
			{ float32_t(      Ex / (1.0 - Vxy * Vyx) * 1e3), float32_t(Ex * Vyx / (1.0 - Vxy * Vyx) * 1e3),          0.0f },
			{ float32_t(Ey * Vxy / (1.0 - Vxy * Vyx) * 1e3), float32_t(      Ey / (1.0 - Vxy * Vyx) * 1e3),          0.0f },
			{                                           0.0,                                          0.0f, float32_t(Es * 1e3) }
		});
	}
};

static constexpr TMaterialProperties MaterialTemplates[] =
{
    // Material parameters from:
    // A Fast Finite Element Solution for Cloth Modelling
    { "Wool"        , 0.26f,  866.0 * 1e-3, 1391.0 * 1e-3, 0.261, 0.162, 0.51 * 1e-3 },
    { "Viscose"     , 0.23f,  245.0 * 1e-3,  366.0 * 1e-3, 0.249, 0.167, 0.38 * 1e-3 },
    { "Polyacrylics", 0.17f, 3057.0 * 1e-3, 1534.0 * 1e-3, 0.150, 0.299, 1.22 * 1e-3 },
    { "Polyester"   , 0.26f, 2400.0 * 1e-3, 3600.0 * 1e-3, 0.250, 0.167, 5.23 * 1e-3 },

    // Material parameters from:
    // Determination of the Elastic Constants of 
    // Plain Woven Fabrics by a Tensile Test in Various Directions
	{ "100% Cotton"        , 0.1503f, 32.559 * 1e3, 12.436 * 1e3, 0.566, 0.216, ShearModulus(0.821, 1.136) * 1e3 },
	{ "100% Wool"          , 0.2348f, 21.860 * 1e3,  8.149 * 1e3, 0.705, 0.263, ShearModulus(0.170, 1.377) * 1e3 },
	{ "95% Wool + 5% lycra", 0.1782f,  0.250 * 1e3,  0.851 * 1e3, 0.071, 0.244, ShearModulus(0.076, 0.599) * 1e3 },
	{ "100% Polyester"     , 0.1646f,  5.152 * 1e3, 11.674 * 1e3, 0.381, 0.864, ShearModulus(0.478, 1.366) * 1e3 },
};

struct TParticle
{
	float32_t mM;
	Eigen::Vector<float32_t, 3> mV;
	Eigen::Vector<float32_t, 3> mXn;
	Eigen::Vector<float32_t, 3> mX0;
    bool mbStationary = false;

	template <bool deformed = true>
	Eigen::Vector<float32_t, 3>& getPosition()
	{
		return deformed ? mXn : mX0;
	}

	template <bool deformed = true>
	const Eigen::Vector<float32_t, 3>& getPosition() const
	{
		return deformed ? mXn : mX0;
	}
};

struct FiniteElementsMethod : public ClothModel
{
	std::vector<TParticle> mParticles;
	std::vector<float32_t> mTrianglesArea;

    float32 mTimeRemainder = 0.0f;

	struct
	{
        float32 timeStep = 0.01f;
		int32 materialTemplate = 0;
        TMaterialProperties materialProperties;
        vec2 rayleighCoefficients = vec2(0.01f, 0.01f);
	} mUserParams;

	FiniteElementsMethod()
	{
        mUserParams.materialProperties = MaterialTemplates[mUserParams.materialTemplate];
		mName = "Finite Elements Method";
	}

	EType getType() const final override
	{
		return EType::FiniteElementsMethod;
	}

	vec3 getVertexPosition(int32 vertexIndex) const final override;

	void init(const vec2 &size, const ivec2 &tessellation) final override;
	void init(const Model *model) final override;
	
	void simulate(ClothSample*, float32_t deltaTime) final override;
	void render(ClothSample*, SampleCallbacks*) final override;

	void onGuiRender(ClothSample*, SampleCallbacks*) final override;
	bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) final override;

private:
	Eigen::Matrix<float32_t, 9, 6> getPlanarProjectedRotation(const trimesh::triangle_t&) const;
	Eigen::Matrix<float32_t, 6, 6> getElementStiffnessMatrix(const Eigen::Vector<float32_t, 2>(&)[3], float32_t) const;
	void assembleCorotatedStiffnessMatrix(Eigen::MatrixXf &K, Eigen::VectorXf &f) const;
	void assembleRayleighDampingMatrix(const Eigen::MatrixXf &M, const Eigen::MatrixXf &K, Eigen::MatrixXf &D) const;
	void assembleParticlesData(Eigen::MatrixXf &M, Eigen::VectorXf &v, Eigen::VectorXf &x) const;
	void assembleExternalForces(Eigen::VectorXf &fext) const;
	float32_t getVoronoiArea(int32 vertexIndex);
	void updateMass();
};

template <bool deformed>
static float32_t getTriangleArea(const trimesh::triangle_t &triangle, const std::vector<TParticle> &particles)
{
	const auto xi = particles[triangle.i].getPosition<deformed>();
	const auto xj = particles[triangle.j].getPosition<deformed>();
	const auto xk = particles[triangle.k].getPosition<deformed>();

	const auto xij = xj - xi;
	const auto xik = xk - xi;

	return xij.cross(xik).norm() / 2.0f;
}

void FiniteElementsMethod::init(const vec2 &size, const ivec2 &tessellation)
{
	std::vector<float3> positions;
	std::vector<float2> texCoords;
	const auto vertexCount = createRectMesh(size, tessellation, positions, texCoords, mTriangles, mMesh);
	
	for (const auto &position : positions)
	{
		TParticle particle;
		particle.mXn[0] = position[0];
		particle.mXn[1] = position[1];
		particle.mXn[2] = position[2];
		particle.mX0 = particle.mXn;
		particle.mV.setZero();
		mParticles.emplace_back(particle);
	}

    mParticles[0].mbStationary = true;
    mParticles[tessellation.x - 1].mbStationary = true;

	for (const auto &triangle : mTriangles)
		mTrianglesArea.emplace_back(getTriangleArea<false>(triangle, mParticles));

	updateMass();
	sharedInit(texCoords);
}

void FiniteElementsMethod::init(const Model *model)
{
	//assert(model->getMeshCount() == 1);
	//assert(model->getMeshInstanceCount(0) == 1);

	//const auto vao = model->getMesh(0)->getVao();

	//assert(vao->getPrimitiveTopology() == Vao::Topology::TriangleList);
	//const auto ib = vao->getIndexBuffer();
	//const auto data = ib->map(Buffer::MapType::Read);
	//for (int32 i = 0; i < ib->getSize() / sizeof(uint8); ++i)
	//{
	//	switch (vao->getIndexBufferFormat())
	//	{
	//	case ResourceFormat::R16Uint:
	//		break;
	//	case ResourceFormat::R32Uint:
	//		EmitTriangle(
	//			reinterpret_cast<uint32*>(data)[i]
	//			reinterpret_cast<uint32*>(data)[i]
	//			reinterpret_cast<uint32*>(data)[i]
	//		);
	//		break;
	//	default:
	//		assert(!"Unsupported index buffer format!");
	//	}
	//}
	//ib->unmap();

	//assert(vao->getVertexBuffersCount() == 1);
	//const auto vb = vao->getVertexBuffer(0);
	//const auto data = vb->map(Buffer::MapType::Read);
	//vb->unmap();
	
	//sharedInit();
}

vec3 FiniteElementsMethod::getVertexPosition(int32 vertexIndex) const
{
	return vec3(
		mParticles[vertexIndex].mXn[0],
		mParticles[vertexIndex].mXn[1],
		mParticles[vertexIndex].mXn[2]
	);
}

void FiniteElementsMethod::updateMass()
{
    const auto &materialProperties = mUserParams.materialProperties;
	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
		mParticles[i].mM = getVoronoiArea(i) * materialProperties.Density;
}

float32_t FiniteElementsMethod::getVoronoiArea(int32 particle)
{
	std::vector<float3> vertices;
	auto faces = mMesh.vertex_face_neighbors(particle);

	auto area = 0.0f;
	for (const auto &face : faces)
	{
		const auto &triangle = mTriangles[face];

		// TODO: most of these these repeat, collect them for all faces
		int32 rotate = -1;
		while (triangle[++rotate] != particle);

		const auto xi = mParticles[triangle[(rotate + 0) % 3]].mX0;
		const auto xj = mParticles[triangle[(rotate + 1) % 3]].mX0;
		const auto xk = mParticles[triangle[(rotate + 2) % 3]].mX0;

		const auto xij = xj - xi;
		const auto xik = xk - xi;
		const auto xjk = xk - xj;

		const auto dir_ij = xij.normalized();
		const auto dir_ik = xik.normalized();
		const auto dir_jk = xjk.normalized();

		const auto cosPhi_i = dir_ij.dot(dir_ik);
		const auto cosPhi_j = -dir_jk.dot(dir_ij);
		const auto cosPhi_k = dir_jk.dot(dir_ik);

		const auto isTriangleObtuse = cosPhi_i < 0.0f || cosPhi_j < 0.0f || cosPhi_k < 0.0f;
		if (!isTriangleObtuse)
		{
			const auto cotgPhi_j = cosPhi_j / sqrtf(1.0f - cosPhi_j * cosPhi_j);
			const auto cotgPhi_k = cosPhi_k / sqrtf(1.0f - cosPhi_k * cosPhi_k);

			area += (xik.dot(xik) * cotgPhi_j + xij.dot(xij) * cotgPhi_k) / 8.0f;
		}
		else
		{
			const auto areaTriangle = mTrianglesArea[face];

			const auto isXiObtuse = cosPhi_i < 0.0f;
			if (isXiObtuse)
				area += areaTriangle / 2.0f;
			else
				area += areaTriangle / 4.0f;
		}
	}
	return area;
}

void FiniteElementsMethod::assembleParticlesData(Eigen::MatrixXf &M, Eigen::VectorXf &v, Eigen::VectorXf &x) const
{
	M.setZero();
	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
	{
		const auto &particle = mParticles[i];

		const auto mass = particle.mM;
		M(3 * i + 0, 3 * i + 0) = mass;
		M(3 * i + 1, 3 * i + 1) = mass;
		M(3 * i + 2, 3 * i + 2) = mass;

		v.segment<3>(3 * i) = particle.mV;
		x.segment<3>(3 * i) = particle.mXn;
	}
}

static Eigen::Matrix<float32_t, 3, 6> getTriangleShapeMatrix(const Eigen::Vector<float32_t, 2>(&triangle)[3], float32_t Ae)
{
	const auto v1 = triangle[0];
	const auto v2 = triangle[1];
	const auto v3 = triangle[2];

	return Eigen::Matrix<float32_t, 3, 6>({
		{          v2.y() - v3.y(),                     0.0f,          v3.y() - v1.y(),                     0.0f,          v1.y() - v2.y(),                     0.0f },
		{                     0.0f,          v3.x() - v2.x(),                     0.0f,          v1.x() - v3.x(),                     0.0f,          v2.x() - v1.x() },
		{ (v3.x() - v2.x()) / 2.0f, (v2.y() - v3.y()) / 2.0f, (v1.x() - v3.x()) / 2.0f, (v3.y() - v1.y()) / 2.0f, (v2.x() - v1.x()) / 2.0f, (v1.y() - v2.y()) / 2.0f }
	}) / (2.0f * Ae);
}

template <bool deformed>
static Eigen::Matrix<float32_t, 2, 3> getPlaneProjectionMatrix(const trimesh::triangle_t &triangle, const std::vector<TParticle> &particles)
{
	const auto x_a = particles[triangle.i].getPosition<deformed>();
	const auto x_b = particles[triangle.j].getPosition<deformed>();
	const auto x_c = particles[triangle.k].getPosition<deformed>();

	const auto n = (x_b - x_a).cross(x_c - x_a);
	const auto px = (x_b - x_a).normalized();
	const auto py = n.cross(px).normalized();

	Eigen::Matrix<float32_t, 2, 3> result;
	result.row(0) = px.transpose();
	result.row(1) = py.transpose();
	
	return result;
}

// 2x2 matrix polar decomposition from:
// Explicit Polar Decomposition and a Near-Characteristic Polynomial: The 2 X 2 Case
Eigen::Matrix<float32_t, 2, 2> getPolarDecompositionRotation(const Eigen::Matrix<float32_t, 2, 2> &M)
{
	const auto &a = M(0, 0);
	const auto &b = M(0, 1);
	const auto &c = M(1, 0);
	const auto &d = M(1, 1);

	if (a * d > b * c)
	{
		const auto beta = sqrtf((a + d) * (a + d) + (b - c) * (b - c));
		return Eigen::Matrix<float32_t, 2, 2>({
			{ a + d, b - c },
			{ c - b, a + d }
		}) / beta;
	}
	else
	{
		const auto beta = sqrtf((a - d) * (a - d) + (b + c) * (b + c));
		return Eigen::Matrix<float32_t, 2, 2>({
			{ a - d, b + c },
			{ b + c, d - a }
		}) / beta;
	}
}

template <bool deformed>
static Eigen::Matrix<float32_t, 2, 2> getCoordinateSystem(const Eigen::Matrix<float32_t, 2, 3> &projection, const trimesh::triangle_t &triangle, const std::vector<TParticle> &particles)
{
	const auto x_a = projection * particles[triangle.i].getPosition<deformed>();
	const auto x_b = projection * particles[triangle.j].getPosition<deformed>();
	const auto x_c = projection * particles[triangle.k].getPosition<deformed>();

	Eigen::Matrix<float32_t, 2, 2> result;
	result.col(0) = x_b - x_a;
	result.col(1) = x_c - x_a;

	return result;
}

Eigen::Matrix<float32_t, 9, 6> FiniteElementsMethod::getPlanarProjectedRotation(const trimesh::triangle_t &triangle) const
{
	const auto P = getPlaneProjectionMatrix<true>(triangle, mParticles);
	const auto P0 = getPlaneProjectionMatrix<false>(triangle, mParticles);

	const auto T = getCoordinateSystem<true>(P, triangle, mParticles);
	const auto S = getCoordinateSystem<false>(P0, triangle, mParticles);

	const auto R = getPolarDecompositionRotation(T * S.inverse());
	const auto PTxR = P.transpose() * R;

	Eigen::Matrix<float32_t, 9, 6> result;
	result.setZero();
	
	result.block<3, 2>(0, 0) = PTxR;
	result.block<3, 2>(3, 2) = PTxR;
	result.block<3, 2>(6, 4) = PTxR;

	return result;
}

Eigen::Matrix<float32_t, 6, 6> FiniteElementsMethod::getElementStiffnessMatrix(const Eigen::Vector<float32_t, 2>(&triangle)[3], float32 A) const
{
    const auto &materialProperties = mUserParams.materialProperties;
	const auto C = materialProperties.getHookeLawMatrix();
	const auto B = getTriangleShapeMatrix(triangle, A);
	return A * B.transpose() * C * B;
}

template <typename TBlockLhs, typename TBlockRhs, typename TMatrix = typename Eigen::Matrix<TBlockRhs::Scalar, TBlockRhs::RowsAtCompileTime, TBlockRhs::ColsAtCompileTime>>
void addBlocks(TBlockLhs &lhs, const TBlockRhs &rhs)
{
	lhs = TMatrix(lhs) + TMatrix(rhs);
}

template <typename TMatrix>
void PrintMatrix(const TMatrix &matrix)
{
	std::string result;
	for (int32 row = 0; row < matrix.rows(); ++row)
	{
		for (int32 col = 0; col < matrix.cols(); ++col)
		{
			result.append(std::to_string(matrix(row, col))).append(" ");
		}
		result.append("\n");
	}
	printToDebugWindow(result);
}

void FiniteElementsMethod::assembleCorotatedStiffnessMatrix(Eigen::MatrixXf &K, Eigen::VectorXf &f0) const
{
	K.setZero();
	f0.setZero();

	for (int32 t = 0, size = int32(mTriangles.size()); t < size; ++t)
	{
		const auto &triangle = mTriangles[t];
		const auto P0 = getPlaneProjectionMatrix<false>(triangle, mParticles);

		Eigen::Vector<float32_t, 2> _x0[3];
		_x0[0] = P0 * mParticles[triangle[0]].mX0;
		_x0[1] = P0 * mParticles[triangle[1]].mX0;
		_x0[2] = P0 * mParticles[triangle[2]].mX0;

		const auto Ke = getElementStiffnessMatrix(_x0, mTrianglesArea[t]);
		const auto R = getPlanarProjectedRotation(triangle);
		const auto RxKe = R * Ke;

		const auto RxKexRT = RxKe * R.transpose();
		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			for (int32 j = 0; j < arraysize(triangle.v); ++j)
			{
				addBlocks(
					K.block<3, 3>(3 * triangle[i], 3 * triangle[j]), RxKexRT.block<3, 3>(3 * i, 3 * j)
				);
			}
		}

		Eigen::Vector<float32_t, 6> x0;
		x0.segment<2>(0) = _x0[0];
		x0.segment<2>(2) = _x0[1];
		x0.segment<2>(4) = _x0[2];
		const auto f0e = -RxKe * x0;

		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			addBlocks(
				f0.segment<3>(3 * triangle[i]), f0e.segment<3>(3 * i)
			);
		}
	}
}

void FiniteElementsMethod::assembleRayleighDampingMatrix(const Eigen::MatrixXf &M, const Eigen::MatrixXf &K, Eigen::MatrixXf &D) const
{
	float32_t alpha = mUserParams.rayleighCoefficients.x, beta = mUserParams.rayleighCoefficients.y;
	D = M * alpha + K * beta;
}

void FiniteElementsMethod::assembleExternalForces(Eigen::VectorXf &fext) const
{
	const auto gravityAcceleration = Eigen::Vector<float32_t, 3>(0.0f, -9.8067f, 0.0f);
	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
	{
		auto &particle = mParticles[i];
		//if (i == 0 || i == 10)
		//{
		//	fext.segment<3>(3 * i) = (mParticles[i].mX0 - mParticles[i].mXn) * 100.0f;
		//}
		//else
		{
			const auto force_gravity = particle.mM * gravityAcceleration;
			fext.segment<3>(3 * i) = force_gravity;
		}
	}
}

void FiniteElementsMethod::simulate(ClothSample*, float32_t delta_t)
{
    //for (mTimeRemainder += delta_t, delta_t = mUserParams.timeStep; mTimeRemainder >= delta_t;)
    { // The model cannot be simulated with interactive framerates.
        //mTimeRemainder -= delta_t;
        if (delta_t == 0.0f) return;
        delta_t = mUserParams.timeStep;
        Eigen::MatrixXf K(3 * mParticles.size(), 3 * mParticles.size());
        Eigen::VectorXf f0(3 * mParticles.size());
        assembleCorotatedStiffnessMatrix(K, f0);

        Eigen::MatrixXf M(3 * mParticles.size(), 3 * mParticles.size());
        Eigen::VectorXf v(3 * mParticles.size());
        Eigen::VectorXf x(3 * mParticles.size());
        assembleParticlesData(M, v, x);

        Eigen::MatrixXf D(3 * mParticles.size(), 3 * mParticles.size());
        assembleRayleighDampingMatrix(M, K, D);

        Eigen::VectorXf fext(3 * mParticles.size());
        assembleExternalForces(fext);

        const Eigen::MatrixXf A = M + delta_t * D + Square(delta_t) * K;
        const Eigen::VectorXf b = -delta_t * (K * x + f0 - fext + delta_t * K * v + D * v);
        const Eigen::VectorXf delta_v = A.colPivHouseholderQr().solve(b); // A * x = b

        for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
        {
            auto &particle = mParticles[i];
            particle.mV += delta_v.segment<3>(3 * i);
            particle.mXn += delta_t * particle.mV;

            if (particle.mbStationary)
            {
                particle.mV.setZero();
                particle.mXn = particle.mX0; // TODO: implement constraints interface
            }
        }
    }
}

void FiniteElementsMethod::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	ClothModel::render(pClothSample, pSample);
}

void FiniteElementsMethod::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

	if (pGui->beginGroup(mName, true))
	{
        pClothSample->addFloatSlider("Const Time Step", mUserParams.timeStep, 0.0001f, 0.1f, false, "%.4f");

		Gui::DropdownList materialTemplateDropdown;
        for (uint32 material = 0u; material < arraysize(MaterialTemplates); ++material)
            materialTemplateDropdown.push_back({ material, MaterialTemplates[material].Name });

        auto &materialProperties = mUserParams.materialProperties;
        if (pGui->addDropdown("Fabric Template", materialTemplateDropdown, reinterpret_cast<uint32&>(mUserParams.materialTemplate)))
        {
            materialProperties = MaterialTemplates[mUserParams.materialTemplate];
            updateMass();
        }

        if (pGui->addFloatSlider("Fabric Density", materialProperties.Density, 0.001f, 10.0f))
            updateMass();
        auto youngModulus = vec2(float32(materialProperties.Ex), float32(materialProperties.Ey));
        if (pGui->addFloat2Slider("Young’s modulus", youngModulus, 0.0f, 100.0f * 1e3f, false))
        {
            materialProperties.Ex = youngModulus.x;
            materialProperties.Ey = youngModulus.y;
        }
        auto poissonRatio = vec2(materialProperties.Vxy, materialProperties.Vyx);
        if (pGui->addFloat2Slider("Poisson’s ratio", poissonRatio, 0.0f, 1.0f, false))
        {
            materialProperties.Vxy = poissonRatio.x;
            materialProperties.Vyx = poissonRatio.y;
        }
        auto shearModulus = float32(materialProperties.Es);
        if (pGui->addFloatSlider("Shear modulus", shearModulus, 0.0f, 1e3f, false))
            materialProperties.Es = shearModulus;

        pGui->addFloat2Slider("Rayleigh Coefficients", mUserParams.rayleighCoefficients, 0.0f, 1.0f, false);

		ClothModel::onGuiRender(pClothSample, pSample);
		pGui->endGroup();
	}
}

bool FiniteElementsMethod::onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&)
{
	return false;
}

template <>
auto ClothModel::createClothModel<ClothModel::FiniteElementsMethod>(const ClothModel *baseParent) -> SharedPtr
{
	auto result = new ::FiniteElementsMethod();
	if (baseParent != nullptr)
	{
		const auto *parent = static_cast<const ::FiniteElementsMethod*>(baseParent);
		result->mUserParams = parent->mUserParams;
	}

	return SharedPtr(result);
}
