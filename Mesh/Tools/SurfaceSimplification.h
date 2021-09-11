#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <queue>
#include <time.h>

#define EPISLON 1E-6F
#define D_PI 6.2831852f

//存储cost，并由Node构建一个小根堆
class Node {
	public:
		int idx;
		double cost;
		Eigen::Vector3d optimized;

		Node(int _idx, double _cost, Eigen::Vector3d _optimized) : idx(_idx), cost(_cost), optimized(_optimized) {};

		//place smallest cost node on the top of heap
		bool operator< (const Node& _rhs) const
		{
			return this->cost > _rhs.cost;
		}
};

namespace DGP {

	class SurfaceSimplification {
	public:
		SurfaceSimplification(MeshKernel::SurfaceMesh& _source)
			: source(_source) 
		{
			this->vertice_count = this->source.VertexSize();
			this->face_count = this->source.FaceSize();
			this->edge_count = this->source.EdgeSize();
		};

		Eigen::Vector3d computeNormal(MeshKernel::Vertex &v0, MeshKernel::Vertex &v1, MeshKernel::Vertex &v2)
		{
			Eigen::Vector3d e10(v1.x() - v0.x(), v1.y() - v0.y(), v1.z() - v0.z());
			Eigen::Vector3d e20(v2.x() - v0.x(), v2.y() - v0.y(), v2.z() - v0.z());
			
			return e20.cross(e10).normalized();
		}

		void computeMatrix_Q();

		void computeCost();

		void decimation();

		void decimation(int n);

		void print_Vector3d(Eigen::Vector3d vec, std::string info)
		{
			std::cout << "info : " + info << " : " << vec.x() << " " << vec.y() << " " << vec.z() << std::endl;
		}

		void print_Matrix4d(Eigen::Matrix4d mat, std::string info)
		{
			std::cout << "Matrix info : " << info << std::endl;
			std::cout << mat.coeffRef(0, 0) << " " << mat.coeffRef(0, 1) << " " << mat.coeffRef(0, 2) << " " << mat.coeffRef(0, 3) << std::endl;
			std::cout << mat.coeffRef(1, 0) << " " << mat.coeffRef(1, 1) << " " << mat.coeffRef(1, 2) << " " << mat.coeffRef(1, 3) << std::endl;
			std::cout << mat.coeffRef(2, 0) << " " << mat.coeffRef(2, 1) << " " << mat.coeffRef(2, 2) << " " << mat.coeffRef(2, 3) << std::endl;
			std::cout << mat.coeffRef(3, 0) << " " << mat.coeffRef(3, 1) << " " << mat.coeffRef(3, 2) << " " << mat.coeffRef(3, 3) << std::endl;

		}

		void test_collapse()
		{
			//this->source.collapse(MeshKernel::EdgeHandle(10));
			//this->source.collapse(MeshKernel::EdgeHandle(20));
			//this->source.collapse(MeshKernel::EdgeHandle(60));
			//this->source.collapse(MeshKernel::EdgeHandle(90));
			//this->source.collapse(MeshKernel::EdgeHandle(100));
			//this->source.collapse(MeshKernel::EdgeHandle(150));
			this->source.collapseTo(MeshKernel::EdgeHandle(257), 1, 1, 1);
		}

		void test_collapse2()
		{
			//this->source.collapse(MeshKernel::EdgeHandle(10));
			//this->source.collapse(MeshKernel::EdgeHandle(20));
			//this->source.collapse(MeshKernel::EdgeHandle(60));
			this->source.collapse2(MeshKernel::EdgeHandle(90));
			//this->source.collapse(MeshKernel::EdgeHandle(100));
			//this->source.collapse(MeshKernel::EdgeHandle(150));
		}

		void test_one_ring_edge()
		{
			MeshKernel::Edge e = this->source.edges(MeshKernel::EdgeHandle(90));
			std::vector<MeshKernel::EdgeHandle> one_ring = this->source.one_ring_edges(e.vh1());
		}

		void test_iteration()
		{
			this->source.test_iteration();
		}

		void test_one_ring(int k)
		{
			std::vector<MeshKernel::EdgeHandle> one_ring_edges;

			MeshKernel::VertexHandle test_vh(0);
			if (k == 1)
			{
				one_ring_edges = this->source.get_one_ring_edges(test_vh);
			}
			else if (k == 2)
			{
				one_ring_edges = this->source.get_one_ring_edges2(test_vh);
			}

			//test get_one_ring_edges results
			std::cout << "start testing the one_ring_edges results" << std::endl;
			int i = 0;
			for (auto eh : one_ring_edges)
			{
				auto e = this->source.edges(eh);
				std::cout << "The " << i++ << " edge : " << e.vh1() << " and " << e.vh2() << std::endl;
			}

			for (auto vh : this->source.NeighborVh(test_vh))
			{
				std::cout << "The neighboring vertexhandle : " << vh.idx() << std::endl;
			}
		}

		void speed_test_one_ring()
		{
			clock_t start, end;
			start = clock();
			for (int i = 0; i < this->source.VertexSize(); i++)
			{
				this->source.get_one_ring_edges(MeshKernel::VertexHandle(i));
			}
			end = clock();
			std::cout << "time comsuption of method 1 : " << end - start << std::endl;

			start = clock();
			for (int i = 0; i < this->source.VertexSize(); i++)
			{
				this->source.get_one_ring_edges2(MeshKernel::VertexHandle(i));
			}
			end = clock();
			std::cout << "time comsuption of method 2 : " << end - start << std::endl;
		}

		

	private:
		MeshKernel::SurfaceMesh& source;
		//MeshKernel::SurfaceMesh& target;
		//MeshKernel::SurfaceMesh& intermediate;
		std::vector<Eigen::Matrix4d> Matrix_Q;
		//std::vector<Eigen::Vector4d> Vector_nd;
		std::priority_queue<Node> queue;
		
		int vertice_count;
		int face_count;
		int edge_count;
	};

}