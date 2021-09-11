#pragma once
#include<iostream>

namespace MeshKernel { 
	class Handle {
	public:
		// 显式构造，不允许 "Handle h = 1;" 的语法进行隐式转换
		explicit Handle(int _idx) : idx_(_idx) {};
		// 赋值构造
		Handle& operator=(const Handle& _h) {
			idx_ = _h.idx_;
			return *this;
		}

		// 基类handle允许用int进行修改
		Handle& operator=(int _idx) {
			idx_ = _idx;
			return *this;
		}

		//判定handle是否存在
		inline bool is_valid() const { return idx_ != -1; }
		/*===========================================handle的比较操作===========================================*/
		inline bool operator<(const Handle& _h) const { return (this->idx_ < _h.idx_); }

		inline bool operator<(int _idx) const { return idx_ < _idx; }

		inline bool operator>(const Handle& _h) const { return (this->idx_ > _h.idx_); }

		inline bool operator>(int _idx) const { return idx_ > _idx; }

		inline bool operator==(const Handle& _h) const { return _h.idx_ == this->idx_; }

		inline bool operator!=(const Handle& _h) const { return _h.idx_ != this->idx_; }

		/*===========================================修改与重置===========================================*/

		inline const int& idx() const { return idx_; }      //取元素

		void idx(const int& _idx) { idx_ = _idx; }      //修改元素

		inline operator int() const { return idx_; } //隐式转换为int

		void reset() { idx_ = -1; }  //初始化

	private:
		int idx_;
	};

	//显式构造能够避免其他类型的handle赋值给当前类型的handle
	class VertexHandle : public Handle { public: explicit VertexHandle(int _idx = -1) : Handle(_idx) {} };
	class EdgeHandle : public Handle { public: explicit EdgeHandle(int _idx = -1) : Handle(_idx) {} };
	class FaceHandle : public Handle { public: explicit FaceHandle(int _idx = -1) : Handle(_idx) {} };
	class CellHandle : public Handle { public: explicit CellHandle(int _idx = -1) : Handle(_idx) {} };
}

/*======================特化各种Handle的哈希映射============================*/
namespace std
{
	template<>
	struct hash<MeshKernel::VertexHandle>
	{
		size_t operator()(const MeshKernel::VertexHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::EdgeHandle>
	{
		size_t operator()(const MeshKernel::EdgeHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::FaceHandle>
	{
		size_t operator()(const MeshKernel::FaceHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::CellHandle>
	{
		size_t operator()(const MeshKernel::CellHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};

}






// note:
// 派生类如果想调用基类已被覆盖的操作符
// 例如 "Handle& operator=(int _idx)"
//void test() {
//	VertexHandle vh(1);
//	vh = 2;                      错误，不允许隐式转换
//	vh.Handle::operator= (2);    正确，调用基类操作符
//}