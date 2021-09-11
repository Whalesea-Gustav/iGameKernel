//#include"IO/IO.h"
//#include"Kernel/Mesh.h"
//#include"Tools/Subdivision.h"
//#include"Tools/Parameterization_CotWeight.h"
//#include"Tools/ARAPShapeInterpolation.h"
//#include <iostream>
//
//#include  <direct.h>  
//#include  <stdio.h> 
//using namespace std;
//
//int main() {
//
//	char   buffer[10000];
//	getcwd(buffer, 10000);
//
//	std::cout << "Working directory is : " << buffer << std::endl;
//
//	MeshKernel::IO io; 
//	//MeshKernel::SurfaceMesh source = io.ReadFile("models/models/example/hw8/0.000000.obj");
//	//MeshKernel::SurfaceMesh target = io.ReadFile("models/models/example/hw8/1.000000.obj");
//
//	MeshKernel::SurfaceMesh source = io.ReadFile("0.000000.obj");
//	MeshKernel::SurfaceMesh target = io.ReadFile("1.000000.obj");
//	MeshKernel::SurfaceMesh intermediate = io.ReadFile("0.000000.obj");
//
//	//DGP::ARAP parameterization(source, target, intermediate);
//	//parameterization.Generate_intermediate_shape(0.3);
//	io.WriteFile(intermediate, "test.obj");
//
//	std::cin.get();
//
//	
//
//	//test1: change vertex position
//	//source.vertices(static_cast<MeshKernel::VertexHandle>(1)).setPosition(10, 10, 10);
//	//io.WriteFile(source, "test.obj");
//
//	//test2: build sparse matrix from triplet
//	
//	//Eigen::SparseMatrix<double> A(2, 2);
//	//std::vector<Eigen::Triplet<double>> triplet;
//	//triplet.push_back(Eigen::Triplet<double>(1, 1, 2));
//	//triplet.push_back(Eigen::Triplet<double>(1, 1, 2));
//	//triplet.push_back(Eigen::Triplet<double>(1, 1, 2));
//	//triplet.push_back(Eigen::Triplet<double>(1, 1, 2));
//	//
//	//A.setFromTriplets(triplet.begin(), triplet.end());
//	//std::cout << std::endl;
//	//std::cout << "triplet test : " << A.coeffRef(1, 1) << std::endl;
//
//	/*MeshKernel::SurfaceMesh newmesh;
//	Subdivision sub(mesh, newmesh);
//	sub.SubTest();
//	io.WriteFile(newmesh,"result1.obj");*/
//
//
//	return 0;
//}