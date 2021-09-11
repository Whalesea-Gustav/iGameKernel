#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#define EPISLON 1E-6F
#define D_PI 6.2831852f

namespace DGP {

	class ARAP {
	public:
		ARAP(MeshKernel::SurfaceMesh& _source, MeshKernel::SurfaceMesh& _target, MeshKernel::SurfaceMesh& _intermediate)
			: source(_source), target(_target), intermediate(_intermediate) {};
		
		void Generate_intermediate_shape(double t);

		void print_Vector3d(Eigen::Vector3d vec, std::string info)
		{
			std::cout << "info : " + info << " : " << vec.x() << " " << vec.y() << " " << vec.z() << std::endl;
		}

	private:
		MeshKernel::SurfaceMesh& source;
		MeshKernel::SurfaceMesh& target;
		MeshKernel::SurfaceMesh& intermediate;
	};
	
}