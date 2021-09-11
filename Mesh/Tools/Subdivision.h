#pragma once
#include"Kernel/Mesh.h"


class Subdivision { 
public:
	Subdivision(MeshKernel::SurfaceMesh& _mesh, MeshKernel::SurfaceMesh& _newmesh) :
	mesh_(_mesh),newmesh_(_newmesh){};
	void SubTest();
private:
	void PosCal();
	void Connect();
	MeshKernel::SurfaceMesh& mesh_;
	MeshKernel::SurfaceMesh& newmesh_;
	std::unordered_map<MeshKernel::EdgeHandle, MeshKernel::Vertex> NewEdgePoint;  // 新的边点
	std::unordered_map<MeshKernel::FaceHandle, MeshKernel::Vertex> NewFacePoint;  // 新的面点
};

void Subdivision::SubTest() {
	PosCal();
	Connect();
}

void Subdivision::PosCal() {
	auto oldmeshvertices = mesh_.allvertices();       // 旧网格上的所有顶点
	for (auto v : oldmeshvertices) {      // 遍历所有旧网格上的顶点
		MeshKernel::Vertex VirtualPoint;             // 当前顶点的虚拟点 根据catmullclark规则得到
		//当前顶点邻接边的另外一个顶点的虚拟点 根据一维三次样条极限点规则得到
		std::unordered_map<MeshKernel::EdgeHandle, MeshKernel::Vertex> D1;   // 检索边得到该虚拟点
		std::unordered_map<MeshKernel::VertexHandle, MeshKernel::Vertex> D2; // 检索点得到该虚拟点
		// 计算该顶点相邻的点的虚拟点
		for (auto nv : mesh_.NeighborVh(v.first)) {        // 遍历当前顶点的所有相邻顶点
			auto eh = mesh_.edgehandle(MeshKernel::Edge(v.first, nv));     // 确定相邻顶点属于哪条相邻边
			std::vector<MeshKernel::VertexHandle> Ci;     // 保存相关点Ci
			std::unordered_set<MeshKernel::VertexHandle> tobeCi;   // 保存Ci的候选者
			// 找到Ci
			for (auto ef : mesh_.NeighborFh(eh)) {          // 遍历该边的两个相邻面
				auto f = mesh_.faces(ef);                   // 当前面
				for (int i = 0; i < f.size(); ++i) {        // 遍历四个顶点
					if (f.vh(i) != v.first) {               // 不是当前处理的点，则加入至Ci候选者
						tobeCi.insert(f.vh(i));             
					}
				}
			}
			for (auto cv : tobeCi) {
				if (mesh_.NeighborVh(nv).count(cv)){      // 找到候选者中符合要求的两个Ci
					Ci.push_back(cv);
				}
			}
			// 计算Di
			MeshKernel::Vertex Di = mesh_.vertices(nv) * 3.0 / 2;   // 3/2 Ei
			for (auto ci : Ci) {
				Di = Di - mesh_.vertices(ci)*1.0 / 4;   // -1/4 Ci
			}
			D1[eh] = MeshKernel::Vertex(Di.x(), Di.y(), Di.z());
			D2[nv] = MeshKernel::Vertex(Di.x(), Di.y(), Di.z());
		}

		int neighborvhsize = mesh_.NeighborEh(v.first).size();      // 得到当前顶点的度数n
		// 计算alpha
		double alpha = (neighborvhsize - 1)*1.0 / (neighborvhsize + 5);
		for (auto f : mesh_.NeighborFh(v.first)) {
			alpha += 4.0 / (neighborvhsize) / (neighborvhsize + 5) * 1.0 / mesh_.faces(f).size();
		}
		// 计算当前顶点的虚拟点
		VirtualPoint = mesh_.vertices(v.first)*1.0 / alpha;         // 第一项
		for (auto fh : mesh_.NeighborFh(v.first)) {                 // 遍历邻接面
			auto f = mesh_.faces(fh);                       
			for (int i = 0; i < f.size(); ++i) {    
				auto curvh = f.vh(i);                              // 当前面的某一顶点
				if (curvh != v.first) {
					auto vpos = mesh_.vertices(curvh);
					if (mesh_.NeighborVh(v.first).count(curvh)) {  //是相邻点则用D替代
						vpos = D2[curvh];
					}
					VirtualPoint = VirtualPoint - vpos* 4.0 / (neighborvhsize)
						/ (neighborvhsize + 5) / alpha / f.size();
				}
				
			}
		}

		// 计算对面点的贡献
		for (auto fh : mesh_.NeighborFh(v.first)) {
			auto f = mesh_.faces(fh);
			MeshKernel::Vertex fi;
			for (int i = 0; i < f.size(); ++i) {
				auto curvh = f.vh(i);                              // 当前面的某一顶点
				auto vpos = mesh_.vertices(curvh);
				if (mesh_.NeighborVh(v.first).count(curvh)) {  //是相邻点则用D替代
					vpos = D2[curvh];
				}
				else if (curvh == v.first) {
					vpos = MeshKernel::Vertex(VirtualPoint.x(), VirtualPoint.y(), VirtualPoint.z());
				}
				fi =  fi + vpos  / f.size();
			}
			if (NewFacePoint.count(fh)) {
				NewFacePoint[fh] = NewFacePoint[fh] + fi;
			}
			else {
				NewFacePoint[fh] = fi;
			}
		}

		// 计算对边点的贡献
		for (auto nv : mesh_.NeighborVh(v.first)) {        // 遍历当前顶点的所有相邻顶点
			auto eh = mesh_.edgehandle(MeshKernel::Edge(v.first, nv));     // 确定相邻顶点属于哪条相邻边
			MeshKernel::Vertex ei;
			ei = VirtualPoint * 1.0 / 3 + D2[nv] * 1.0 / 6;
			for (auto ef : mesh_.NeighborFh(eh)) {          // 遍历该边的两个相邻面
				auto f = mesh_.faces(ef);                   // 当前面
				for (int i = 0; i < f.size(); ++i) {        // 遍历所有顶点
					if (f.vh(i) != v.first &&f.vh(i) != nv) {
						if (mesh_.NeighborVh(nv).count(f.vh(i))) {
							ei = ei + mesh_.vertices(f.vh(i))*1.0 / 12;
						}
						if (mesh_.NeighborVh(v.first).count(f.vh(i))) {
							ei = ei + mesh_.vertices(f.vh(i))*1.0 / 12;
						}
					}
				}
			}
			if (NewEdgePoint.count(eh)) {
				NewEdgePoint[eh] = NewEdgePoint[eh] + ei;
			}
			else {
				NewEdgePoint[eh] = ei;
			}

		}

	}
}

void Subdivision::Connect() {
	auto oldmeshallface = mesh_.allfaces();
	for (auto f : oldmeshallface) {
		auto facepoint = NewFacePoint[f.first] *1.0/ f.second.size();
		auto facepointhandle = newmesh_.AddVertex(MeshKernel::Vertex(facepoint.x(), facepoint.y(), facepoint.z()));

		for (int i = 0; i <4; ++i) {
			auto vh = f.second.vh(i);
			auto oldv = mesh_.vertices(vh);
			auto newpointhandle = newmesh_.AddVertex(MeshKernel::Vertex(oldv.x(), oldv.y(), oldv.z()));
			MeshKernel::Vertex edgepoint1 = NewEdgePoint[f.second.eh((i + 1) % f.second.size())] / 2.0;
			MeshKernel::Vertex edgepoint2 = NewEdgePoint[f.second.eh(i)]/2.0;
			
			auto edgepointhandle1 = newmesh_.AddVertex(MeshKernel::Vertex(edgepoint1.x(), edgepoint1.y(), edgepoint1.z()));
			auto edgepointhandle2 = newmesh_.AddVertex(MeshKernel::Vertex(edgepoint2.x(), edgepoint2.y(), edgepoint2.z()));
			newmesh_.AddFace(std::vector<MeshKernel::VertexHandle>{newpointhandle, edgepointhandle1, facepointhandle, edgepointhandle2});
			
		}
	}
}