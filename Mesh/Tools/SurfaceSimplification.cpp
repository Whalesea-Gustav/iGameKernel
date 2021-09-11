#include "SurfaceSimplification.h"
#include <queue>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
//#include <queue>
#include <iostream>

using namespace MeshKernel;
void DGP::SurfaceSimplification::computeMatrix_Q()
{
	Matrix_Q.clear();

	for (int i = 0; i < this->vertice_count; i++)
	{
		Matrix_Q.push_back(Eigen::Matrix4d::Zero());
	}

	//遍历所有的面

	for (auto fp : this->source.allfaces())
	{
		FaceHandle fh = fp.first;
		//std::cout << fh.idx() << std::endl;
		Face f = fp.second;
		//calculate triangle Q
		VertexHandle vh0 = f.vh(0);
		VertexHandle vh1 = f.vh(1);
		VertexHandle vh2 = f.vh(2);

		Vertex &v0 = this->source.vertices(vh0);
		Vertex &v1 = this->source.vertices(vh1);
		Vertex &v2 = this->source.vertices(vh2);

		Eigen::Vector3d ev0(v0.x(), v0.y(), v0.z());
		Eigen::Vector3d ev1(v1.x(), v1.y(), v1.z());
		Eigen::Vector3d ev2(v2.x(), v2.y(), v2.z());

		Eigen::Vector3d normal = computeNormal(v0, v1, v2);
		//std::cout << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;

		Eigen::Vector4d n_bar(normal.x(), normal.y(), normal.z(), 0);
		n_bar.coeffRef(3) = -1.0 * normal.dot(ev0);
		Eigen::Matrix4d Qi = n_bar * n_bar.transpose();
		
		Matrix_Q[static_cast<int>(vh0)] += Qi;
		Matrix_Q[static_cast<int>(vh1)] += Qi;
		Matrix_Q[static_cast<int>(vh2)] += Qi;
	}
	//test
	/*
	for (int i = 0; i < this->source.VertexSize(); i++)
	{
		print_Matrix4d(Matrix_Q[i], std::to_string(i));
	}
	*/
}

void DGP::SurfaceSimplification::computeCost()
{
	//std::priority_queue<Node> heap;
	queue = std::priority_queue<Node>();

	for (auto ep : this->source.alledges())
	{
		EdgeHandle eh = ep.first;
		Edge &e = ep.second;

		VertexHandle &vh1 = e.vh1();
		VertexHandle &vh2 = e.vh2();

		Eigen::Matrix4d Q_Edge = Matrix_Q[static_cast<int>(vh1)] + Matrix_Q[static_cast<int>(vh2)];
	
		// Compute the optimal contraction target 

		//Eigen::Matrix3d A = Q_Edge.block<3, 3>(0, 0);
		Eigen::Matrix3d A;
		A << Q_Edge.coeffRef(0, 0), Q_Edge.coeffRef(0, 1), Q_Edge.coeffRef(0, 2),
			Q_Edge.coeffRef(1, 0), Q_Edge.coeffRef(1, 1), Q_Edge.coeffRef(1, 2),
			Q_Edge.coeffRef(2, 0), Q_Edge.coeffRef(2, 1), Q_Edge.coeffRef(2, 2);
		//Eigen::Vector3d b = Q_Edge.block<3, 1>(3, 0);
		Eigen::Vector3d b;
		b(0) = Q_Edge.coeffRef(3, 0);
		b(1) = Q_Edge.coeffRef(3, 1);
		b(2) = Q_Edge.coeffRef(3, 2);

		Eigen::Vector3d x = A.inverse() * (-1 * b);
		//std::cout << "info : optimized point : " << x.x() << " " << x.y() << " " << x.z() << std::endl;

		Vertex v1 = this->source.vertices(vh1);
		Vertex v2 = this->source.vertices(vh2);
		
		Eigen::Vector3d ev1(v1.x(), v1.y(), v1.z());
		Eigen::Vector3d ev2(v2.x(), v2.y(), v2.z());

		//test
		//if (eh.idx() == 3)
		//{
		//	print_Vector3d(ev1, "vertex 1");
		//	print_Vector3d(ev2, "vertex 2");
		//	print_Vector3d(x, "optimized point");
		//}

		Eigen::Vector4d x_bar(x.x(), x.y(), x.z(), 1);
		double cost = x_bar.transpose() * Q_Edge * x_bar;
		//std::cout << "cost : " << cost << std::endl;

		Node node(eh.idx(), cost, x);

		queue.push(node);

	}
}

void DGP::SurfaceSimplification::decimation()
{
	//visualize the queue
	//for (auto ep : this->source.alledges())
	//{
	//	Node largest = queue.top();
	//	queue.pop();
	//	std::cout << "idx : " << largest.idx << " |  cost : " << largest.cost << std::endl;
	//}

	//to do : 
	//Iteratively remove the pair(v1, v2) of least cost from the heap, 
	//contract this pair, 
	//and update the costs of all valid pairs involving v1.

	Node node= queue.top();
	queue.pop();
	//std::cout << "top idx : " << node.idx << std::endl;
	this->source.collapseTo(MeshKernel::EdgeHandle(node.idx), node.optimized.x(), node.optimized.y(), node.optimized.z());
	computeMatrix_Q();
	computeCost();

}

void DGP::SurfaceSimplification::decimation(int n)
{
	for (int i = 0; i < n; i++)
	{
		decimation();
	}
}
