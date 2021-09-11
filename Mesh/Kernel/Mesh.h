#pragma once
#include"Cell.h"

namespace MeshKernel { 
	// 网格的概念
	class Mesh {
	public:
		/*=========================读写元素===============================*/
			// 读取ID为i的顶点
		Vertex& vertices(VertexHandle _vh);
		const Vertex vertices(VertexHandle _vh) const;        // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用   
		// 读取ID为i的边
		Edge& edges(EdgeHandle _eh);
		const Edge& edges(EdgeHandle _eh) const;
		// 读取ID为i的面
		Face& faces(FaceHandle _fh);
		const Face faces(FaceHandle _fh) const;

		const std::unordered_map<VertexHandle, Vertex>& allvertices() const { return vertices_; }
		const std::unordered_map<EdgeHandle, Edge>& alledges() const { return edges_; }
		const std::unordered_map<FaceHandle, Face>& allfaces() const { return faces_; }

		/*====================根据元素得到对应ID=========================*/
		const VertexHandle vertexhanle(Vertex _vertex) const;
		const EdgeHandle edgehandle(Edge& _edge) const;
		const FaceHandle facehandle(Face& _face) const;
		/*======================得到邻接关系============================*/
		// 顶点的邻接点
		std::unordered_set<VertexHandle> NeighborVh(VertexHandle _vh);
		// 顶点的邻接边
		std::unordered_set<EdgeHandle> NeighborEh(VertexHandle _vh);
		// 顶点的邻接面
		std::unordered_set<FaceHandle> NeighborFh(VertexHandle _vh);
		// 边的邻接边
		std::unordered_set<EdgeHandle> NeighborEh(EdgeHandle _eh);
		// 边的邻接面
		std::unordered_set<FaceHandle> NeighborFh(EdgeHandle _eh);
		// 面的邻接面
		std::unordered_set<FaceHandle> NeighborFh(FaceHandle _fh);
		/*=========================添加元素=============================*/
		VertexHandle AddVertex(const Vertex& _v);
		EdgeHandle AddEdge(const VertexHandle& _vh1, const VertexHandle& _vh2);
		virtual FaceHandle AddFace(const std::vector<VertexHandle>& _vhs);

		/*=========================删除元素=============================*/
		// 删除低级元素时删除一切邻接的高级元素
		// 删除高级元素时保留仍被其他元素使用的低级元素
		VertexHandle DeleteVertex(const VertexHandle& _vh);
		EdgeHandle DeleteEdge(const EdgeHandle& _eh);
		FaceHandle DeleteFace(const FaceHandle& _fh);



		/*=========================生成唯一ID========================*/
		VertexHandle GenVertexHandle() { return (VertexHandle)VertexHandleID_++; }
		EdgeHandle GenEdgeHandle() { return (EdgeHandle)EdgeHandleID_++; }
		FaceHandle GenFaceHandle() { return (FaceHandle)FaceHandleID_++; }

		size_t VertexSize() { return vertices_.size(); }
		size_t FaceSize() { return faces_.size(); }
		size_t EdgeSize() { return edges_.size(); }
	protected:
		/*============下一个可使用的ID=========*/
		int VertexHandleID_ = 0;
		int EdgeHandleID_ = 0;
		int FaceHandleID_ = 0;

		/*=============修改邻接关系============*/
		// 将该面添加至面所包含的点和边的相邻面中
		void AddFace2Neighbor(const FaceHandle& _f);
		// 将该边添加至边所包含的点相邻边中
		void AddEdge2Neighbor(const EdgeHandle& _e);
		// Todo: Delete Neighbor
		void DeleteFace2Neighbor(const FaceHandle& _fh);
		void DeleteEdge2Neighbor(const EdgeHandle& _eh);

	protected:
		// handle到元素的对应
		std::unordered_map<VertexHandle, Vertex> vertices_;
		std::unordered_map<EdgeHandle, Edge> edges_;
		std::unordered_map<FaceHandle, Face> faces_;

		// 元素到handle的对应
		std::unordered_map<Vertex, VertexHandle> Vertex2Vh_;
		std::unordered_map<Edge, EdgeHandle> Edge2Eh_;
		std::unordered_map<Face, FaceHandle> Face2Fh_;

		// 邻接关系
		std::unordered_map<VertexHandle, std::unordered_set<EdgeHandle>> NeighborEhOfVertex_;          //点的邻接边
		std::unordered_map<VertexHandle, std::unordered_set<FaceHandle>> NeighborFhOfVertex_;          //点的邻接面
		std::unordered_map<EdgeHandle, std::unordered_set<FaceHandle>> NeighborFhOfEdge_;              //边的邻接面
	protected:
		virtual void InitMesh(const std::vector<Vertex>& _vertices,
			const std::vector<std::vector<VertexHandle>>& _elements) = 0;                              // 必须重写
		Mesh& operator=(const Mesh& _mesh);
	};

	// 曲面网格
	class SurfaceMesh : public Mesh {
	public:
		/*==========================构造函数==============================*/
		SurfaceMesh() {};
		SurfaceMesh(const std::vector<Vertex>& _vertices, const std::vector<std::vector<VertexHandle>>& _faces) {
			InitMesh(_vertices, _faces);
		}
		//SurfaceMesh(const SurfaceMesh& _surfacemesh);
		/*=============初始化网格=============*/
		virtual void InitMesh(const std::vector<Vertex>& _vertices,
			const std::vector<std::vector<VertexHandle>>& _elements) override;
		SurfaceMesh& operator=(const SurfaceMesh& _surfacemesh);

		/*=========================检验是否在边界上=============================*/
		bool IsOnBoundary(VertexHandle vh);
		bool SurfaceMesh::IsOnBoundary(EdgeHandle vh);

		/*=========================折叠翻转算法=============================*/
		void collapse(EdgeHandle& _eh);
		void collapseTo(EdgeHandle& _eh, double x, double y, double z);
		void collapse2(EdgeHandle& _eh);
		void flip(EdgeHandle& _eh);

		/*=========================1-邻域对边=============================*/
		void test_iteration();
		EdgeHandle& findEdge(VertexHandle& from_v, VertexHandle& to_v);
		std::vector<EdgeHandle> one_ring_edges(VertexHandle& _vh); //查找1-领域所有的点构成的环结构
		std::vector<EdgeHandle> get_one_ring_edges(VertexHandle& vh);
		std::vector<EdgeHandle> get_one_ring_edges2(VertexHandle& vh);
	};

	class VolumeMesh : public Mesh {
	};

	class TriMesh : SurfaceMesh{};
	class QuadMesh : SurfaceMesh{};

	class TetMesh : VolumeMesh{};
	class HexMesh : VolumeMesh{};
}