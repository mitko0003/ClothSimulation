#include "Eigen/Eigen"
#include "External/HalfEdge/trimesh.h"

#include "ClothSample.h"

static constexpr float32_t ShearModulus(float32_t E45, float32_t V45)
{
	return E45 * 0.5f / (1.0f + V45);
}

struct TElasticProperties // For orthotropic materials
{
	float32_t Ex, Ey;    // Young’s modulus in MPa
	float32_t Vxy, Vyx;  // Poisson’s ratio
	float32_t Es;        // Shear modulus

	Eigen::Matrix<float32_t, 3, 3> getHookeLawMatrix() const
	{
		return Eigen::Matrix<float32_t, 3, 3>({
			{       Ex / (1.0f - Vxy * Vyx), Ex * Vyx / (1.0f - Vxy * Vyx), 0.0f },
			{ Ey * Vxy / (1.0f - Vxy * Vyx),       Ey / (1.0f - Vxy * Vyx), 0.0f },
			{                          0.0f,                          0.0f,   Es }
		});
	}
};

// Experimental values from:
// Determination of the Elastic Constants of 
// Plain Woven Fabrics by a Tensile Test in Various Directions
static constexpr TElasticProperties ElasticProperties[] =
{
	{ 32.559f, 12.436f, 0.566f, 0.243f, ShearModulus(0.821f, 1.136f) }, // 100% Cotton
	{ 21.860f,  8.149f, 0.705f, 0.277f, ShearModulus(0.170f, 1.377f) }, // 100% Wool
	{  0.250f,  0.851f, 0.071f, 0.196f, ShearModulus(0.076f, 0.599f) }, // 95% Wool + 5% lycra
	{  5.152f, 11.674f, 0.381f, 0.779f, ShearModulus(0.478f, 1.366f) }, // 100% Polyester
};

struct TParticle
{
	float32_t mM;
	Eigen::Vector<float32_t, 3> mV;
	Eigen::Vector<float32_t, 3> mXn;
	Eigen::Vector<float32_t, 3> mX0;

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

	TElasticProperties mElasticProperties;

	FiniteElementsMethod()
	{
		mName = "Finite Elements Method";
	}

	EType getType() const final override
	{
		return EType::FiniteElementsMethod;
	}

	vec3 getVertexPosition(int32 vertexIndex) const final override;
	void init(int32 width, int32 height) final override;
	void init(const Model *model) final override;
	
	void simulate(float32_t deltaTime) final override;
	void render(ClothSample*, SampleCallbacks*) final override;

	void onGuiRender(ClothSample*, SampleCallbacks*) final override;
	bool onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&) final override;

private:
	Eigen::Matrix<float32_t, 9, 6> getPlanarProjectedRotation(const trimesh::triangle_t&) const;
	Eigen::Matrix<float32_t, 6, 6> getElementStiffnessMatrix(const trimesh::triangle_t&, float32_t) const;
	void assembleCorotatedStiffnessMatrix(Eigen::MatrixXf &K, Eigen::VectorXf &f) const;
	void assembleRayleighDampingMatrix(const Eigen::MatrixXf &M, const Eigen::MatrixXf &K, Eigen::MatrixXf &D) const;
	void assembleParticlesData(Eigen::MatrixXf &M, Eigen::VectorXf &v, Eigen::VectorXf &x) const;
	void assembleExternalForces(Eigen::VectorXf &fext) const;
	float32_t getVoronoiArea(int32 vertexIndex);
};

template <bool deformed>
static float32_t getTriangleArea(const trimesh::triangle_t &triangle, const std::vector<TParticle> &particles)
{
	const auto xi = particles[triangle.i].getPosition<deformed>();
	const auto xj = particles[triangle.j].getPosition<deformed>();
	const auto xk = particles[triangle.k].getPosition<deformed>();

	const auto xij = xj - xi;
	const auto xik = xk - xi;

	return xij.cross(xik)[0] / 2.0f;
}

void FiniteElementsMethod::init(int32 width, int32 height)
{
	std::vector<float3> positions;
	const auto vertexCount = createRectMesh(width, height, positions, mTriangles, mMesh);
	
	for (const auto &position : positions)
	{
		TParticle particle;
		particle.mXn[0] = position[0];
		particle.mXn[1] = position[1];
		particle.mXn[2] = position[2];
		particle.mX0 = particle.mXn;
		particle.mV = particle.mV.Zero();
		mParticles.emplace_back(particle);
	}

	for (const auto &triangle : mTriangles)
		mTrianglesArea.emplace_back(getTriangleArea<false>(triangle, mParticles));

	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
		mParticles[i].mM = getVoronoiArea(i) * mMaterialDensity;

	sharedInit();
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

float32_t FiniteElementsMethod::getVoronoiArea(int32 particle)
{
	std::vector<float3> vertices;
	auto faces = mMesh.vertex_face_neighbors(particle);

	auto area = 0.0f;
	for (const auto &face : faces) // TODO: this returns indexes in mTriangles
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

template <bool deformed>
static Eigen::Matrix<float32_t, 3, 6> getTriangleShapeMatrix(const trimesh::triangle_t &triangle, const std::vector<TParticle> &particles)
{
	const auto v1 = particles[triangle.i].getPosition<deformed>();
	const auto v2 = particles[triangle.j].getPosition<deformed>();
	const auto v3 = particles[triangle.k].getPosition<deformed>();

	return Eigen::Matrix<float32_t, 3, 6>({
		{          v2.y() - v3.y(),                     0.0f,          v3.y() - v1.y(),                     0.0f,          v1.y() - v2.y(),                     0.0f },
		{                     0.0f,          v3.x() - v2.x(),                     0.0f,          v1.x() - v3.x(),                     0.0f,          v2.x() - v1.x() },
		{ (v3.x() - v2.x()) / 2.0f, (v2.y() - v3.y()) / 2.0f, (v1.x() - v3.x()) / 2.0f, (v3.y() - v1.y()) / 2.0f, (v2.x() - v1.x()) / 2.0f, (v1.y() - v2.y()) / 2.0f }
	});
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
	result = result.Zero();
	
	result.block<3, 2>(0, 0) = PTxR;
	result.block<3, 2>(3, 2) = PTxR;
	result.block<3, 2>(6, 4) = PTxR;

	return result;
}

Eigen::Matrix<float32_t, 6, 6> FiniteElementsMethod::getElementStiffnessMatrix(const trimesh::triangle_t &triangle, float32 A) const
{
	const auto C = mElasticProperties.getHookeLawMatrix();
	const auto B = getTriangleShapeMatrix<true>(triangle, mParticles);
	return A * B.transpose() * C * B;
}

template <typename TBlockLhs, typename TBlockRhs, typename TMatrix = typename Eigen::Matrix<TBlockRhs::Scalar, TBlockRhs::RowsAtCompileTime, TBlockRhs::ColsAtCompileTime>>
void addBlocks(TBlockLhs &lhs, const TBlockRhs &rhs)
{
	lhs = TMatrix(lhs) + TMatrix(rhs);
}

void FiniteElementsMethod::assembleCorotatedStiffnessMatrix(Eigen::MatrixXf &K, Eigen::VectorXf &f) const
{
	for (int32 t = 0, size = int32(mTriangles.size()); t < size; ++t)
	{
		const auto &triangle = mTriangles[t];

		const auto Ke = getElementStiffnessMatrix(triangle, mTrianglesArea[t]);
		const auto R = getPlanarProjectedRotation(triangle);
		const auto RxKe = R * Ke;

		const auto RxKexRT = RxKe * R.transpose();
		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			for (int32 j = 0; j < arraysize(triangle.v); ++j)
			{
				addBlocks(
					K.block<3, 3>(3 * triangle.v[i], 3 * triangle.v[j]), RxKexRT.block<3, 3>(3 * i, 3 * j)
				);
			}
		}

		const auto P0 = getPlaneProjectionMatrix<false>(triangle, mParticles);
		
		Eigen::Vector<float32_t, 6> x0;
		x0.segment<2>(0) = P0 * mParticles[triangle[0]].mX0;
		x0.segment<2>(2) = P0 * mParticles[triangle[1]].mX0;
		x0.segment<2>(4) = P0 * mParticles[triangle[2]].mX0;
		const auto f0 = RxKe * x0;

		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			addBlocks(
				f.segment<3>(3 * triangle.v[i]), f0.segment<3>(3 * i)
			);
		}
	}
}

void FiniteElementsMethod::assembleRayleighDampingMatrix(const Eigen::MatrixXf &M, const Eigen::MatrixXf &K, Eigen::MatrixXf &D) const
{
	float32_t alpha = 0.5f, beta = 0.5f;
	D = M * alpha + K * beta;
}

void FiniteElementsMethod::assembleExternalForces(Eigen::VectorXf &fext) const
{
	const auto gravityAcceleration = Eigen::Vector<float32_t, 3>(0.0f, -9.8067f, 0.0f);
	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
	{
		auto &particle = mParticles[i];
		const auto force_gravity = particle.mM * gravityAcceleration;
		fext.segment<3>(3 * i) = force_gravity;
	}
}

void FiniteElementsMethod::simulate(float32_t delta_t)
{
	Eigen::MatrixXf K(3 * mParticles.size(), 3 * mParticles.size());
	Eigen::VectorXf f(3 * mParticles.size());
	assembleCorotatedStiffnessMatrix(K, f);

	Eigen::MatrixXf M(3 * mParticles.size(), 3 * mParticles.size());
	Eigen::VectorXf v(3 * mParticles.size());
	Eigen::VectorXf x(3 * mParticles.size());
	assembleParticlesData(M, v, x);

	Eigen::MatrixXf D(3 * mParticles.size(), 3 * mParticles.size());
	assembleRayleighDampingMatrix(M, K, D);

	Eigen::VectorXf fext(3 * mParticles.size());
	assembleExternalForces(fext);

	const Eigen::MatrixXf A = M + delta_t * D + Square(delta_t) * K;
	const Eigen::VectorXf b = K * x + f - fext + delta_t * K * v + D * v;
	const Eigen::VectorXf delta_v = A.colPivHouseholderQr().solve(b); // A * x = b

	for (int32 i = 0, size = int32(mParticles.size()); i < size; ++i)
	{
		auto &particle = mParticles[i];
		particle.mV += delta_v.segment<3>(3 * i);
		particle.mXn += delta_t * particle.mV;
	}

	mParticles[0].mXn = mParticles[0].mX0; // TODO: implement constraints interface
}

void FiniteElementsMethod::render(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	ClothModel::render(pClothSample, pSample);
}

void FiniteElementsMethod::onGuiRender(ClothSample *pClothSample, SampleCallbacks *pSample)
{
	auto *pGui = pSample->getGui();

	if (pGui->beginGroup(mName))
	{
		ClothModel::onGuiRender(pClothSample, pSample);
		pGui->endGroup();
	}
}

bool FiniteElementsMethod::onMouseEvent(ClothSample*, SampleCallbacks*, const MouseEvent&)
{
	return false;
}

template <>
ClothModel* ClothModel::createClothModel<ClothModel::FiniteElementsMethod>()
{
	return new ::FiniteElementsMethod();
}