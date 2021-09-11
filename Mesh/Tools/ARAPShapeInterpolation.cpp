#include "ARAPShapeInterpolation.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace DGP {

	void ARAP::Generate_intermediate_shape(double t) {

		//1.svd分解并存储Affine transformation矩阵（Jacobian矩阵）
		size_t face_count = this->source.FaceSize();
		size_t vertex_count = this->source.VertexSize();
		std::cout << "vertex numbers : " << vertex_count << std::endl;

		std::vector<Eigen::Vector3d> xx(face_count);
		std::vector<Eigen::Vector3d> yy(face_count);
		std::vector<double> area(face_count);
		std::vector<Eigen::Matrix2d> S(face_count);
		std::vector<double> angle(face_count);

		Eigen::Matrix2d I = Eigen::MatrixXd::Identity(2, 2);

		std::vector<std::vector<MeshKernel::VertexHandle>> v_id(face_count);

		for (auto f_pair : this->source.allfaces())
		{
			MeshKernel::Face f = f_pair.second;
			MeshKernel::FaceHandle fh = f_pair.first;
			//int id = fh.idx();
			MeshKernel::VertexHandle vh0 = f.vh(0);
			MeshKernel::VertexHandle vh1 = f.vh(1);
			MeshKernel::VertexHandle vh2 = f.vh(2);
			
			v_id[fh.idx()].push_back(vh0);
			v_id[fh.idx()].push_back(vh1);
			v_id[fh.idx()].push_back(vh2);

			Eigen::Vector3d e2t0(this->source.vertices(vh2).x() - this->source.vertices(vh0).x(), this->source.vertices(vh2).y() - this->source.vertices(vh0).y(), 0);
			Eigen::Vector3d e1t0(this->source.vertices(vh1).x() - this->source.vertices(vh0).x(), this->source.vertices(vh1).y() - this->source.vertices(vh0).y(), 0);
			//面积 也用于求逆矩阵
			area[fh.idx()] = 0.5 * e2t0.cross(e1t0).norm();

			//xx, yy用于求逆矩阵
			xx[fh.idx()] = Eigen::Vector3d(this->source.vertices(vh2).x() - this->source.vertices(vh1).x(), 
				this->source.vertices(vh0).x() - this->source.vertices(vh2).x(),
				this->source.vertices(vh1).x() - this->source.vertices(vh0).x());

			yy[fh.idx()] = Eigen::Vector3d(this->source.vertices(vh1).y() - this->source.vertices(vh2).y(),
				this->source.vertices(vh2).y() - this->source.vertices(vh0).y(),
				this->source.vertices(vh0).y() - this->source.vertices(vh1).y());

			//target网格的坐标，用于求逆矩阵
			Eigen::Vector3d u(this->target.vertices(vh0).x(), this->target.vertices(vh1).x(), this->target.vertices(vh2).x());
			Eigen::Vector3d v(this->target.vertices(vh0).y(), this->target.vertices(vh1).y(), this->target.vertices(vh2).y());

			
			Eigen::Matrix2d J;
			J << yy[fh.idx()].dot(u) / (2 * area[fh.idx()]), xx[fh.idx()].dot(u) / (2 * area[fh.idx()]),
				yy[fh.idx()].dot(v) / (2 * area[fh.idx()]), xx[fh.idx()].dot(v) / (2 * area[fh.idx()]);
			//或者直接求逆来计算

			//svd decomposition
			Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);

			//define scaling
			Eigen::Matrix2d sigma;
			sigma << svd.singularValues()[0], 0.0,
				0.0, svd.singularValues()[1];
			S[fh.idx()] = svd.matrixV() * sigma * svd.matrixV().transpose();
			//std::cout << "Scaling Matrix : " << std::endl;
			//std::cout << S[fh.idx()](0, 0) << S[fh.idx()](0, 1) << std::endl;
			//std::cout << S[fh.idx()](1, 0) << S[fh.idx()](1, 1) << std::endl;
			//define rotation: 2D rotation matrix
			//[cos, -sin]
			//[sin,  cos]
			Eigen::Matrix2d R = svd.matrixU() * svd.matrixV().transpose();
			angle[fh.idx()] = atan2(R.coeffRef(1, 0), R.coeffRef(1, 1));
			//std::cout << "angle : " << angle[fh.idx()] << std::endl;
		}

		//fixing point
		double u0 = this->source.vertices(static_cast<MeshKernel::VertexHandle>(vertex_count - 1)).x();
		double v0 = this->source.vertices(static_cast<MeshKernel::VertexHandle>(vertex_count - 1)).y();


		//using t to interpolate Jacobian matrix

		Eigen::SparseMatrix<double> A(2 * vertex_count - 2, 2 * vertex_count - 2);

		std::vector<Eigen::Triplet<double>> triplet;
		
		for (auto f_pair : this->source.allfaces())
		{
			auto fh = f_pair.first;
			int id = fh.idx();
			//auto f = f_pair.second;

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if ((v_id[id][i] == vertex_count - 1) || (v_id[id][j] == vertex_count - 1))
					{
						continue;
					}
					else
					{
						triplet.push_back(Eigen::Triplet<double>(2 * v_id[id][i], 2 * v_id[id][j], yy[id][i] * yy[id][j] / (area[id] * area[id] * 2)));
						triplet.push_back(Eigen::Triplet<double>(2 * v_id[id][i], 2 * v_id[id][j], xx[id][i] * xx[id][j] / (area[id] * area[id] * 2)));
						triplet.push_back(Eigen::Triplet<double>(2 * v_id[id][i] + 1, 2 * v_id[id][j] + 1, yy[id][i] * yy[id][j] / (area[id] * area[id] * 2)));
						triplet.push_back(Eigen::Triplet<double>(2 * v_id[id][i] + 1, 2 * v_id[id][j] + 1, xx[id][i] * xx[id][j] / (area[id] * area[id] * 2)));
					}

				}
			}	
		}
		A.setFromTriplets(triplet.begin(), triplet.end());
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		if (solver.info() == Eigen::Success)
		{
			std::cout << "solver factorize success" << std::endl;
		}
		else
		{
			std::cout << "solver factorize failed" << std::endl;
			std::cout << solver.lastErrorMessage() << std::endl;
		}

		Eigen::VectorXd b(2 * vertex_count - 2);
		b.setZero();
		for (auto f_pair : this->source.allfaces())
		{
			auto fh = f_pair.first;
			int id = fh.idx();
			//auto f = f_pair.second;

			Eigen::Matrix2d M, R;
			double angle_t = angle[id] * t;
			//插值得到t时刻的放射变换矩阵
			R << cos(angle_t), -sin(angle_t),
				sin(angle_t), cos(angle_t);
			M = R * (I * (1 - t) + S[id] * t);

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if ((v_id[id][i] == vertex_count - 1) || (v_id[id][j] == vertex_count - 1))
					{
						if ((v_id[id][i] != vertex_count - 1) && (v_id[id][j] == vertex_count - 1))
						{
							b[2 * v_id[id][i]] += -yy[id][i] * yy[id][j] * u0 / (area[id] * area[id] * 2);
							b[2 * v_id[id][i]] += -xx[id][i] * xx[id][j] * u0 / (area[id] * area[id] * 2);
							b[2 * v_id[id][i] + 1] += -yy[id][i] * yy[id][j] * v0 / (area[id] * area[id] * 2);
							b[2 * v_id[id][i] + 1] += -xx[id][i] * xx[id][j] * v0 / (area[id] * area[id] * 2);
						}
					}
				}

				if (v_id[id][i] != vertex_count - 1)
				{
					b[2 * v_id[id][i]] += yy[id][i] * M(0, 0) / area[id];
					b[2 * v_id[id][i]] += xx[id][i] * M(0, 1) / area[id];
					b[2 * v_id[id][i] + 1] += yy[id][i] * M(1, 0) / area[id];
					b[2 * v_id[id][i] + 1] += xx[id][i] * M(1, 1) / area[id];
				}
			}
		}
		std::cout << "building b success" << std::endl;

		Eigen::VectorXd result = solver.solve(b);

		std::cout << "solving AX=b success" << std::endl;

		for (int i = 0; i < vertex_count - 1; i++)
		{
			this->intermediate.vertices(static_cast<MeshKernel::VertexHandle>(i)).setPosition(result[2 * i], result[2 * i + 1], 0);
		}
	}
}