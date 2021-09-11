#pragma once
#include"Face.h"

namespace MeshKernel { 
	class Cell {
	public:
		/*=========================构造函数=============================*/
		Cell() :n_(0) {}
		Cell(const std::vector<VertexHandle>& _vertices, const std::vector<EdgeHandle>& _edges,
			 const std::vector<FaceHandle>& _faces){
			assert(_faces.size() + _vertices.size() - _edges.size() == 2);     //是否满足欧拉公式
			assert(_faces.size() >= 4 && _vertices.size() >= 4);
			n_ = _faces.size();
			vertices_.assign(_vertices.begin(), _vertices.end());
			edges_.assign(_edges.begin(), _edges.end());
			faces_.assign(_faces.begin(), _faces.end());
		}
		Cell(const Cell& _c) {
			*this = _c;
		}
		/*=======================基本操作符==========================*/
		inline Cell& operator=(const Cell& _c) {
			n_ = _c.n_;
			vertices_.assign(_c.vertices_.begin(), _c.vertices_.end());
			edges_.assign(_c.edges_.begin(), _c.edges_.end());
			faces_.assign(_c.faces_.begin(), _c.faces_.end());
			return *this;
		}
		bool operator==(const Cell& _c) const {
			return (n_ == _c.n_)  && isSameFaces(_c);
		}
		/*=======================读写元素==========================*/
		const VertexHandle vh(int k) const { return vertices_[k]; }          // 常量对象只支持读
		const EdgeHandle eh(int k) const { return edges_[k]; }
		const FaceHandle fh(int k) const { return faces_[k]; }
		
		VertexHandle& vh(int k) { return vertices_[k]; }                     // 非常量对象支持读写
		EdgeHandle& eh(int k) { return edges_[k]; }
		FaceHandle& fh(int k) { return faces_[k]; }
		
		size_t size() const { return n_; }
		/*================得到有序的顶点handle排列==========================*/
		std::vector<VertexHandle> getSortedVertexHandle() const {
			std::vector<VertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}
	private:
		// 判断输入体的面是否与当前体全部相同
		bool isSameFaces(const Cell& _c) const {
			std::unordered_set<FaceHandle> fhset;
			for (const auto& f1 : _c.faces_) fhset.insert(f1);
			for (const auto& f2 : _c.faces_) fhset.insert(f2);
			return fhset.size() == n_;
		}
	private:
		// 保存该体所有点的handle
		std::vector<VertexHandle> vertices_;
		// 保存该体所有边的handle
		std::vector<EdgeHandle> edges_;
		// 保存该体所有面的handle
		std::vector<FaceHandle> faces_;
		// 保存该体的面数
		size_t n_;                     
	};

}

/*======================特化Cell的哈希映射============================*/
// To do: 修改Hash映射
namespace std
{
	template<> struct hash<MeshKernel::Cell>
	{
		size_t operator()(const MeshKernel::Cell& _c)const
		{
			size_t res = 0;
			assert(_c.size() >= 4);
			auto cv = _c.getSortedVertexHandle();
			for (int i = 0; i < cv.size(); ++i) {
				res ^= hash<int>()(cv[i]);
			}
			return res;
		}
	};
}