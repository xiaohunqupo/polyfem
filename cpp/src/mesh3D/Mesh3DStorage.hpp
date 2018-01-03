#ifndef MESH_STORAGE_HPP__
#define MESH_STORAGE_HPP__

#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

namespace poly_fem
{
	struct Vertex
	{
		int id;
		std::vector<double> v;
		std::vector<uint32_t> neighbor_vs;
		std::vector<uint32_t> neighbor_es;
		std::vector<uint32_t> neighbor_fs;
		std::vector<uint32_t> neighbor_hs;

		bool boundary;
		bool boundary_hex;
	};
	struct Edge
	{
		int id;
		std::vector<uint32_t> vs;
		std::vector<uint32_t> neighbor_fs;
		std::vector<uint32_t> neighbor_hs;

		bool boundary;
		bool boundary_hex;
	};
	struct Face
	{
		int id;
		std::vector<uint32_t> vs;
		std::vector<uint32_t> es;
		std::vector<uint32_t> neighbor_hs;
		bool boundary;
		bool boundary_hex;
	};

	struct Element
	{
		int id;
		std::vector<uint32_t> vs;
		std::vector<uint32_t> es;
		std::vector<uint32_t> fs;
		std::vector<bool> fs_flag;
		bool hex = false;
		std::vector<double> v_in_Kernel;
	};

	enum MeshType 
	{
		Tri = 0,
		Qua,
		HSur,
		Tet,
		Hyb,
		Hex
	};

	struct Mesh3DStorage
	{
		MeshType type;
		Eigen::MatrixXd points;
		std::vector<Vertex> vertices;
		std::vector<Edge> edges;
		std::vector<Face> faces;
		std::vector<Element> elements;
	};

	struct Mesh_Quality
	{
		std::string Name;
		double min_Jacobian;
		double ave_Jacobian;
		double deviation_Jacobian;
		VectorXd V_Js;
		VectorXd H_Js;
		VectorXd Num_Js;
		int32_t V_num, H_num;
	};
}

#endif
