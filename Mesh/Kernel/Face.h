#pragma once
#include"Edge.h"
#include<algorithm>
namespace MeshKernel{
	// 网格中的面片
	class Face {
	public: 
		/*=========================构造函数=============================*/
		Face() :n_(0) {}
		Face(const std::vector<VertexHandle>& _vertices, const std::vector<EdgeHandle>& _edges)  {
			assert(_edges.size() == _vertices.size());
			assert(_edges.size() >= 3);
			n_ = _vertices.size();
			edges_.assign(_edges.begin(), _edges.end());
			vertices_.assign(_vertices.begin(), _vertices.end());
		}
		Face(const Face& _f) {
			*this = _f;
		}
		/*=======================基本操作符==========================*/
		inline Face& operator=(const Face& _f) {
			n_ = _f.n_;
			edges_.assign(_f.edges_.begin(), _f.edges_.end());
			vertices_.assign(_f.vertices_.begin(), _f.vertices_.end());
			return *this;
		}
		bool operator==(const Face& _f) const {
			return (n_== _f.n_) && isSameEdges(_f);
		}

		/*=======================读写元素==========================*/
		const VertexHandle vh(int k) const { return vertices_[k]; }            // 常量对象只支持读
		const EdgeHandle eh(int k) const { return edges_[k]; }
		VertexHandle& vh(int k) { return vertices_[k]; }		               // 非常量对象支持读写
		EdgeHandle& eh(int k) { return edges_[k]; }
		
		const size_t size() const { return n_; }

		/*================得到有序的顶点handle排列==========================*/
		std::vector<VertexHandle> getSortedVertexHandle() const {
			std::vector<VertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}
	private:
		// 判断输入面的边是否与当前面全部相同
		bool isSameEdges(const Face& _f) const {
			std::unordered_set<EdgeHandle> ehset;
			for (const auto& e1 : _f.edges_) ehset.insert(e1);
			for (const auto& e2 : edges_) ehset.insert(e2);
			return ehset.size() == n_;
		}
	private:
		// 保存该面所有边的handle
		std::vector<EdgeHandle> edges_;
		// 保存该面所有点的handle
		std::vector<VertexHandle> vertices_;
		// 保存该面的边数
		size_t n_;
	};	
}

/*======================特化Face的哈希映射============================*/
// To do: 修改Hash映射
namespace std
{
	template<> struct hash <MeshKernel::Face>
	{
		size_t operator()(const MeshKernel::Face& _f)const
		{
			size_t res = 0;
			assert(_f.size() >= 3);
			auto fv = _f.getSortedVertexHandle();
			for (int i = 0; i < fv.size(); ++i) {
				res ^= hash<int>()(fv[i]);
			}
			return res;
		}
	};
}