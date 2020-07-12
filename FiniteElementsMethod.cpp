#include "FiniteElementsMethod.h"

static constexpr float ShearModulus(float E45, float V45)
{
	return E45 * 0.5f / (1.0f + V45);
}

struct TElasticProperties // For orthotropic materials
{
	float Ex, Ey;    // Young’s modulus in MPa
	float Vxy, Vyx;  // Poisson’s ratio
	float Es;        // Shear modulus

	Eigen::Matrix<float, 3, 3> getHookeLawMatrix() const
	{
		return Eigen::Matrix<float, 3, 3>(
			      Ex / (1.0f - Vxy * Vyx), Ex * Vyx / (1.0f - Vxy * Vyx), 0.0f,
			Ey * Vxy / (1.0f - Vxy * Vyx),       Ey / (1.0f - Vxy * Vyx), 0.0f,
		                             0.0f,                          0.0f,   Es
		);
	}
};

// Experimental values from:
// Determination of the Elastic Constants of 
// Plain Woven Fabrics by a Tensile Test in Various Directions
static constexpr TElasticProperties ElasticProperties[] =
{
	{ 32.559f, 12.436f, 0.566f, 0.243, ShearModulus(0.821f, 1.136f) }, // 100% Cotton
	{ 21.860f,  8.149f, 0.705f, 0.277, ShearModulus(0.170f, 1.377f) }, // 100% Wool
	{  0.250f,  0.851f, 0.071f, 0.196, ShearModulus(0.076f, 0.599f) }, // 95% Wool + 5% lycra
	{  5.152f, 11.674f, 0.381f, 0.779, ShearModulus(0.478f, 1.366f) }, // 100% Polyester
};

struct FiniteElementsMethod : public ClothModel
{
	std::vector<Eigen::Vector<float, 3>> mXn;
	std::vector<Eigen::Vector<float, 3>> mX0;

	std::vector<trimesh::triangle_t> mTriangles;
	trimesh::trimesh_t mMesh;

	TElasticProperties mElasticProperties;

	Buffer::SharedPtr mpVBPositions;
	Buffer::SharedPtr mpVBNormals;
	Buffer::SharedPtr mpVBColors;
	Buffer::SharedPtr mpIndexBuffer;

	FiniteElementsMethod()
	{
		mName = "Finite Elements Method";
	}

	void init(int32_t width, int32_t height) final override;
	void init(const Model *model) final override;
	void simulate(float deltaTime) final override;
	void render() final override;

private:
	Eigen::Matrix<float, 9, 6> getPlanarProjectedRotation(const trimesh::triangle_t&) const;
	Eigen::Matrix<float, 6, 6> getElementStiffnessMatrix(const trimesh::triangle_t&) const;
	void assembleCorotatedStiffnessMatrix() const;
	void updateVoronoiArea(int32 vertexIndex);
};

void FiniteElementsMethod::updateVoronoiArea(int32 vertexIndex)
{
	std::vector<float3> vertices;
	auto faces = mMesh.vertex_face_neighbors(vertexIndex);

	auto area = 0.0f;
	for (int32 face = 0; face < faces.size() / 3; ++face) // TODO: this returns indexes in mTriangles
	{
		// TODO: most of these these repeat, collect them for all faces
		const auto xi = vertices[faces[face + 0]];
		const auto xj = vertices[faces[face + 1]];
		const auto xk = vertices[faces[face + 2]];

		const auto xij = xj - xi;
		const auto xik = xk - xi;
		const auto xjk = xk - xj;

		const auto dir_ij = normalize(xij);
		const auto dir_ik = normalize(xik);
		const auto dir_jk = normalize(xjk);

		const auto cosPhi_i = dot(dir_ij, dir_ik);
		const auto cosPhi_j = dot(dir_jk, dir_ij);
		const auto cosPhi_k = dot(dir_jk, dir_ik);

		const auto isTriangleObtuse = cosPhi_i < 0.0f || cosPhi_j < 0.0f || cosPhi_k < 0.0f;
		if (!isTriangleObtuse)
		{
			const auto cotgPhi_j = cosPhi_j / sqrtf(1.0f - cosPhi_j * cosPhi_j);
			const auto cotgPhi_k = cosPhi_k / sqrtf(1.0f - cosPhi_k * cosPhi_k);

			area += (dot(xik, xik) * cotgPhi_j + dot(xjk, xij) * cotgPhi_k) / 8.0f;
		}
		else
		{
			auto areaTriangle = 0.0f;

			const auto isXiObtuse = cosPhi_i < 0.0f;
			if (isXiObtuse)
				area += areaTriangle / 2.0f;
			else
				area += areaTriangle / 4.0f;
		}
	}
}

Eigen::Matrix<float, 3, 6> getTriangleShapeMatrix(const trimesh::triangle_t &triangle, const std::vector<Eigen::Vector<float, 3>> &vertices)
{
	const auto v1 = vertices[triangle.i];
	const auto v2 = vertices[triangle.j];
	const auto v3 = vertices[triangle.k];

	return Eigen::Matrix<float, 3, 6>(
		         v2.y - v3.y,                 0.0f,          v3.y - v1.y,                 0.0f,          v1.y - v2.y,                 0.0f,
		                0.0f,          v3.x - v2.x,                 0.0f,          v1.x - v3.x,                 0.0f,          v2.x - v1.x,
		(v3.x - v2.x) / 2.0f, (v2.y - v3.y) / 2.0f, (v1.x - v3.x) / 2.0f, (v3.y - v1.y) / 2.0f, (v2.x - v1.x) / 2.0f, (v1.y - v2.y) / 2.0f
	);
}

float getTriangleArea(const trimesh::triangle_t &triangle, const std::vector<Eigen::Vector<float, 3>> &vertices)
{
	return 0.0f;
}

Eigen::Matrix<float, 2, 3> getPlaneProjectionMatrix(const trimesh::triangle_t &triangle, const std::vector<Eigen::Vector<float, 3>> &vertices)
{
	const auto x_a = vertices[triangle.i];
	const auto x_b = vertices[triangle.j];
	const auto x_c = vertices[triangle.k];

	const auto n = (x_b - x_a).cross(x_c - x_a);
	const auto px = (x_b - x_a).normalize();
	const auto py = n.cross(px).normalize();

	Eigen::Matrix<float, 2, 3> result;
	result.row(0) = px.transpose();
	result.row(1) = py.transpose();
	
	return result;
}

// 2x2 matrix polar decomposition from:
// Explicit Polar Decomposition and a Near-Characteristic Polynomial: The 2 X 2 Case
Eigen::Matrix<float, 2, 2> getPolarDecompositionRotation(const Eigen::Matrix<float, 2, 2> &M)
{
	const auto &a = M(0, 0);
	const auto &b = M(0, 1);
	const auto &c = M(1, 0);
	const auto &d = M(1, 1);

	if (a * d > b * c)
	{
		const auto beta = sqrtf((a + d) * (a + d) + (b - c) * (b - c));
		return Eigen::Matrix<float, 2, 2>(
			a + d, b - c,
			c - b, a + d
		) / beta;
	}
	else
	{
		const auto beta = sqrtf((a - d) * (a - d) + (b + c) * (b + c));
		return Eigen::Matrix<float, 2, 2>(
			a - d, b + c,
			b + c, d - a
		) / beta;
	}
}

Eigen::Matrix<float, 2, 2> getCoordinateSystem(const Eigen::Matrix<float, 2, 3> &projection, const trimesh::triangle_t &triangle, const std::vector<Eigen::Vector<float, 3>> &vertices)
{
	const auto x_a = projection * vertices[triangle.i];
	const auto x_b = projection * vertices[triangle.j];
	const auto x_c = projection * vertices[triangle.k];

	Eigen::Matrix<float, 2, 2> result;
	result.col(0) = x_b - x_a;
	result.col(1) = x_c - x_a;

	return result;
}

Eigen::Matrix<float, 9, 6> FiniteElementsMethod::getPlanarProjectedRotation(const trimesh::triangle_t &triangle) const
{
	const auto P = getPlaneProjectionMatrix(triangle, mXn);
	const auto P0 = getPlaneProjectionMatrix(triangle, mX0);

	const auto T = getCoordinateSystem(P, triangle, mXn);
	const auto S = getCoordinateSystem(P0, triangle, mX0);

	const auto R = getPolarDecompositionRotation(T * S.inverse());
	const auto RxPT = R * P.transpose();

	auto result = Eigen::Matrix<float, 9, 6>::Zero();
	result.block<3, 3>(0, 0) = RxPT;
	result.block<3, 3>(3, 3) = RxPT;
	result.block<3, 3>(6, 6) = RxPT;

	return result;
}

Eigen::Matrix<float, 6, 6> FiniteElementsMethod::getElementStiffnessMatrix(const trimesh::triangle_t &triangle) const
{
	const auto C = mElasticProperties.getHookeLawMatrix();
	const auto B = getTriangleShapeMatrix(triangle, mXn);
	const auto A = getTriangleArea(triangle, mX0);
	return A * B.transpose() * C * B;
}

void FiniteElementsMethod::assembleCorotatedStiffnessMatrix() const
{
	Eigen::MatrixXf K(3 * mXn.size(), 3 * mXn.size());
	Eigen::VectorXf f(3 * mXn.size());

	for (const auto &triangle : mTriangles)
	{
		const auto K = getElementStiffnessMatrix(triangle);
		const auto R = getPlanarProjectedRotation(triangle);
		const auto RxK = R * K;

		const auto RxKxRT = RxK * R.transpose();
		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			for (int32 j = 0; j < arraysize(triangle.v); ++j)
			{
				RxKxRT.block<3, 3>(3 * i, 3 * j).addTo(
					K.block<3, 3>(3 * triangle.v[i], 3 * triangle.v[j])
				);
			}
		}

		const auto P0 = getPlaneProjectionMatrix(triangle, mX0);
		
		Eigen::Vector<float, 6> x0;
		x0.block<2, 1>(0, 0) = P0 * mX0[triangle[0]];
		x0.block<2, 1>(2, 0) = P0 * mX0[triangle[1]];
		x0.block<2, 1>(4, 0) = P0 * mX0[triangle[2]];
		const auto f0 = RxK * x0;

		for (int32 i = 0; i < arraysize(triangle.v); ++i)
		{
			f0.block<3, 1>(3 * i, 1).addTo(
				f.block<3, 1>(3 * triangle.v[i], 1)
			);
		}
	}
}

void FiniteElementsMethod::init(const Model *model)
{	
	assert(model->getMeshCount() == 1);
	assert(model->getMeshInstanceCount(0) == 1);

	const auto vao = model->getMesh(0)->getVao();

	assert(vao->getPrimitiveTopology() == Vao::Topology::TriangleList);
	const auto ib = vao->getIndexBuffer();
	const auto data = ib->map(Buffer::MapType::Read);
	for (int32 i = 0; i < ib->getSize() / sizeof(uint8); ++i)
	{
		switch (vao->getIndexBufferFormat())
		{
		case ResourceFormat::R16Uint:
			break;
		case ResourceFormat::R32Uint:
			EmitTriangle(
				reinterpret_cast<uint32*>(data)[i]
				reinterpret_cast<uint32*>(data)[i]
				reinterpret_cast<uint32*>(data)[i]
			);
			break;
		default:
			assert(!"Unsupported index buffer format!");
		}
	}
	ib->unmap();

	assert(vao->getVertexBuffersCount() == 1);
	const auto vb = vao->getVertexBuffer(0);
	const auto data = vb->map(Buffer::MapType::Read);
	vb->unmap();
}

void FiniteElementsMethod::simulate(float deltaTime)
{

}

void FiniteElementsMethod::render()
{

}

FiniteElementsMethod* getFiniteElementsMethodCloth()
{
	return new FiniteElementsMethod();
}