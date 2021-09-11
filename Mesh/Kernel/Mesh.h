#pragma once
#include"Cell.h"

namespace MeshKernel { 
	// ����ĸ���
	class Mesh {
	public:
		/*=========================��дԪ��===============================*/
			// ��ȡIDΪi�Ķ���
		Vertex& vertices(VertexHandle _vh);
		const Vertex vertices(VertexHandle _vh) const;        // unordered_map �� [] ���������ǳ�����Ա�������޷��Գ�������ʹ��   
		// ��ȡIDΪi�ı�
		Edge& edges(EdgeHandle _eh);
		const Edge& edges(EdgeHandle _eh) const;
		// ��ȡIDΪi����
		Face& faces(FaceHandle _fh);
		const Face faces(FaceHandle _fh) const;

		const std::unordered_map<VertexHandle, Vertex>& allvertices() const { return vertices_; }
		const std::unordered_map<EdgeHandle, Edge>& alledges() const { return edges_; }
		const std::unordered_map<FaceHandle, Face>& allfaces() const { return faces_; }

		/*====================����Ԫ�صõ���ӦID=========================*/
		const VertexHandle vertexhanle(Vertex _vertex) const;
		const EdgeHandle edgehandle(Edge& _edge) const;
		const FaceHandle facehandle(Face& _face) const;
		/*======================�õ��ڽӹ�ϵ============================*/
		// ������ڽӵ�
		std::unordered_set<VertexHandle> NeighborVh(VertexHandle _vh);
		// ������ڽӱ�
		std::unordered_set<EdgeHandle> NeighborEh(VertexHandle _vh);
		// ������ڽ���
		std::unordered_set<FaceHandle> NeighborFh(VertexHandle _vh);
		// �ߵ��ڽӱ�
		std::unordered_set<EdgeHandle> NeighborEh(EdgeHandle _eh);
		// �ߵ��ڽ���
		std::unordered_set<FaceHandle> NeighborFh(EdgeHandle _eh);
		// ����ڽ���
		std::unordered_set<FaceHandle> NeighborFh(FaceHandle _fh);
		/*=========================���Ԫ��=============================*/
		VertexHandle AddVertex(const Vertex& _v);
		EdgeHandle AddEdge(const VertexHandle& _vh1, const VertexHandle& _vh2);
		virtual FaceHandle AddFace(const std::vector<VertexHandle>& _vhs);

		/*=========================ɾ��Ԫ��=============================*/
		// ɾ���ͼ�Ԫ��ʱɾ��һ���ڽӵĸ߼�Ԫ��
		// ɾ���߼�Ԫ��ʱ�����Ա�����Ԫ��ʹ�õĵͼ�Ԫ��
		VertexHandle DeleteVertex(const VertexHandle& _vh);
		EdgeHandle DeleteEdge(const EdgeHandle& _eh);
		FaceHandle DeleteFace(const FaceHandle& _fh);



		/*=========================����ΨһID========================*/
		VertexHandle GenVertexHandle() { return (VertexHandle)VertexHandleID_++; }
		EdgeHandle GenEdgeHandle() { return (EdgeHandle)EdgeHandleID_++; }
		FaceHandle GenFaceHandle() { return (FaceHandle)FaceHandleID_++; }

		size_t VertexSize() { return vertices_.size(); }
		size_t FaceSize() { return faces_.size(); }
		size_t EdgeSize() { return edges_.size(); }
	protected:
		/*============��һ����ʹ�õ�ID=========*/
		int VertexHandleID_ = 0;
		int EdgeHandleID_ = 0;
		int FaceHandleID_ = 0;

		/*=============�޸��ڽӹ�ϵ============*/
		// ��������������������ĵ�ͱߵ���������
		void AddFace2Neighbor(const FaceHandle& _f);
		// ���ñ���������������ĵ����ڱ���
		void AddEdge2Neighbor(const EdgeHandle& _e);
		// Todo: Delete Neighbor
		void DeleteFace2Neighbor(const FaceHandle& _fh);
		void DeleteEdge2Neighbor(const EdgeHandle& _eh);

	protected:
		// handle��Ԫ�صĶ�Ӧ
		std::unordered_map<VertexHandle, Vertex> vertices_;
		std::unordered_map<EdgeHandle, Edge> edges_;
		std::unordered_map<FaceHandle, Face> faces_;

		// Ԫ�ص�handle�Ķ�Ӧ
		std::unordered_map<Vertex, VertexHandle> Vertex2Vh_;
		std::unordered_map<Edge, EdgeHandle> Edge2Eh_;
		std::unordered_map<Face, FaceHandle> Face2Fh_;

		// �ڽӹ�ϵ
		std::unordered_map<VertexHandle, std::unordered_set<EdgeHandle>> NeighborEhOfVertex_;          //����ڽӱ�
		std::unordered_map<VertexHandle, std::unordered_set<FaceHandle>> NeighborFhOfVertex_;          //����ڽ���
		std::unordered_map<EdgeHandle, std::unordered_set<FaceHandle>> NeighborFhOfEdge_;              //�ߵ��ڽ���
	protected:
		virtual void InitMesh(const std::vector<Vertex>& _vertices,
			const std::vector<std::vector<VertexHandle>>& _elements) = 0;                              // ������д
		Mesh& operator=(const Mesh& _mesh);
	};

	// ��������
	class SurfaceMesh : public Mesh {
	public:
		/*==========================���캯��==============================*/
		SurfaceMesh() {};
		SurfaceMesh(const std::vector<Vertex>& _vertices, const std::vector<std::vector<VertexHandle>>& _faces) {
			InitMesh(_vertices, _faces);
		}
		//SurfaceMesh(const SurfaceMesh& _surfacemesh);
		/*=============��ʼ������=============*/
		virtual void InitMesh(const std::vector<Vertex>& _vertices,
			const std::vector<std::vector<VertexHandle>>& _elements) override;
		SurfaceMesh& operator=(const SurfaceMesh& _surfacemesh);

		/*=========================�����Ƿ��ڱ߽���=============================*/
		bool IsOnBoundary(VertexHandle vh);
		bool SurfaceMesh::IsOnBoundary(EdgeHandle vh);

		/*=========================�۵���ת�㷨=============================*/
		void collapse(EdgeHandle& _eh);
		void collapseTo(EdgeHandle& _eh, double x, double y, double z);
		void collapse2(EdgeHandle& _eh);
		void flip(EdgeHandle& _eh);

		/*=========================1-����Ա�=============================*/
		void test_iteration();
		EdgeHandle& findEdge(VertexHandle& from_v, VertexHandle& to_v);
		std::vector<EdgeHandle> one_ring_edges(VertexHandle& _vh); //����1-�������еĵ㹹�ɵĻ��ṹ
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