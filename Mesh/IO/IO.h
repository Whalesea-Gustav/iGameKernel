#pragma
#include"Kernel/Mesh.h"
#include<fstream>
#include<string>
#include<sstream>


namespace MeshKernel {
	class IO {
	public:
		SurfaceMesh ReadFile(const std::string& _InputFile);
		bool WriteFile(const SurfaceMesh& _mesh,const std::string& _OutputFile);

	private:
		void ReOrderVertexHandle(const SurfaceMesh& _mesh);
		std::vector<VertexHandle> reorderedvh_;                        // 重排顶点
		std::unordered_map<VertexHandle, std::size_t> newvh_;          // 新的顶点handle
	};
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadFile(const std::string& _InputFile) {
	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<Vertex> vertices;
	std::vector<std::vector<VertexHandle>> faces;
	std::string line;

	std::cout << "Reading " << _InputFile << " File" << std::endl;
	// std::cout << inputfile.good() << std::endl;
	int i = 0;
	while (inputfile) {
		line.clear();
		getline(inputfile, line);
		if (line[0] == '#') {
			continue;// 如果是注释的情况
		}
		std::stringstream linestream;
		linestream.str(line);

		char flag = 0;
		linestream >> flag;
		if (flag == 'v') {
			double x, y, z;
			linestream >> x >> y >> z;
			vertices.push_back(Vertex(x, y, z));
		}
		else if (flag == 'f') {
			// 首先确定该面的边数
			int n = 0;
			for (int i = 1; i < line.size(); i++) if (line[i] == ' ') n++;
			// 创建面
			// std::cout << " 加入的面的点数 : " << n << std::endl;
			std::vector<VertexHandle> face(n);
			for (int i = 0; i < n; i++) {
				int vh;
				linestream >> vh;// 加入每一个点的序号
				face[i] = (VertexHandle)(vh - 1);
			}
			faces.push_back(face);
		}
	}
	return SurfaceMesh(vertices, faces);
	inputfile.close();
}

bool MeshKernel::IO::WriteFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	std::ofstream outputfile(_OutputFile, std::ios::out);
	std::cout << "Writing " << _OutputFile << " File" << std::endl;
	ReOrderVertexHandle(_mesh);
	for (VertexHandle vh : reorderedvh_) {
		Vertex v(_mesh.vertices(vh));
		outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		outputfile << "f";
		for (int i = 0; i < f.second.size(); ++i) {
			outputfile << " " << newvh_[f.second.vh(i)]+1;
		}
		outputfile << std::endl;
	}
	outputfile.close();
	return true;
}

void MeshKernel::IO::ReOrderVertexHandle(const SurfaceMesh& _mesh) {
	auto allv = _mesh.allvertices();
	int idx = 0;
	for (auto v : allv) {
		reorderedvh_.push_back(v.first);
		newvh_[v.first] = idx++;
	}
}