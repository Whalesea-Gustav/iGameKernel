#pragma once
#include"Handle.h"
#include"Vertex.h"
#include<iostream>


namespace MeshKernel { 
	// 网格中的边
	class Edge {
	public:
		/*=========================构造函数=============================*/
		Edge() :vertices_(2, (VertexHandle)-1) {}
		Edge(const Edge& _e) {
			vertices_ = _e.vertices_;
		}
		Edge(const VertexHandle& _vh1, const VertexHandle& _vh2) :
			vertices_{ _vh1,_vh2 } {
		}

		/*=======================基本操作符==========================*/
		inline Edge& operator=(const Edge& _e) {
			vertices_ = _e.vertices_;
			return *this;
		}
		inline bool operator==(const Edge& _e) const {
			return ((vh1() == _e.vh1() && vh2() == _e.vh2())
				|| (vh1() == _e.vh2() && vh2() == _e.vh1()));
		}

		/*=========================读写元素==============================*/
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

		/*=================对当前vertices_进行排序=======================*/
		inline void SortVertices() {
			if (vertices_[0] > vertices_[1]) {
				std::swap(vertices_[0], vertices_[1]);
			}
		}

	private:
		// 保存该边两个顶点的vertexhandle
		std::vector<VertexHandle> vertices_;


	};

}

/*======================特化Edge的哈希映射============================*/
// To do: 修改Hash映射
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