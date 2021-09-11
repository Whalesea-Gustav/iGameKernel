#pragma once
#include"Handle.h"
#include"Vertex.h"
#include<iostream>


namespace MeshKernel { 
	// �����еı�
	class Edge {
	public:
		/*=========================���캯��=============================*/
		Edge() :vertices_(2, (VertexHandle)-1) {}
		Edge(const Edge& _e) {
			vertices_ = _e.vertices_;
		}
		Edge(const VertexHandle& _vh1, const VertexHandle& _vh2) :
			vertices_{ _vh1,_vh2 } {
		}

		/*=======================����������==========================*/
		inline Edge& operator=(const Edge& _e) {
			vertices_ = _e.vertices_;
			return *this;
		}
		inline bool operator==(const Edge& _e) const {
			return ((vh1() == _e.vh1() && vh2() == _e.vh2())
				|| (vh1() == _e.vh2() && vh2() == _e.vh1()));
		}

		/*=========================��дԪ��==============================*/
		const VertexHandle vh1() const { return vertices_[0]; }
		const VertexHandle vh2() const { return vertices_[1]; }
		VertexHandle& vh1() { return vertices_[0]; }
		VertexHandle& vh2() { return vertices_[1]; }

		const VertexHandle vh(int k) const { 
			assert(k <= 1 && k >= 0);
			return vertices_[k]; 
		}
		VertexHandle& vh(int k) {
			assert(k <= 1 && k >= 0);
			return vertices_[k];
		}

		/*=================�Ե�ǰvertices_��������=======================*/
		inline void SortVertices() {
			if (vertices_[0] > vertices_[1]) {
				std::swap(vertices_[0], vertices_[1]);
			}
		}

	private:
		// ����ñ����������vertexhandle
		std::vector<VertexHandle> vertices_;


	};

}

/*======================�ػ�Edge�Ĺ�ϣӳ��============================*/
// To do: �޸�Hashӳ��
namespace std
{
	template<> struct hash<MeshKernel::Edge>
	{
		size_t operator()(const MeshKernel::Edge& _e)const
		{
			if (_e.vh1() <= _e.vh2()) {
				return hash<int>()(_e.vh1()) ^ hash<int>()(_e.vh2());
			}
			else {
				return hash<int>()(_e.vh2()) ^ hash<int>()(_e.vh1());
			}
			
		}
	};
}