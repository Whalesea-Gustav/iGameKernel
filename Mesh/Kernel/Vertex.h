#pragma once
#include<vector>
#include<iostream>
#include<vector>
#include<unordered_set>
#include<unordered_map>
#include <cassert>
#include "GeometryKernel.h"
#include "Handle.h"

namespace MeshKernel{  
	// 网格中的顶点
	class Vertex {
	public:
		Vertex() :Vec3(3, 0.0) {}
		Vertex(const Vertex& _v) :Vec3(_v.Vec3) {}
		Vertex(double x, double y, double z) :Vec3{ x,y,z } {}

		/*===========读写单一分量元素=================*/
		// 常量成员只能读元素
		const double x() const { return Vec3[0]; }
		const double y() const { return Vec3[1]; }
		const double z() const { return Vec3[2]; }

		// 普通成员可以读写元素
		double& x() { return Vec3[0]; }
		double& y() { return Vec3[1]; }
		double& z() { return Vec3[2]; }

		/*============基本运算=================*/
		inline Vertex operator+(const Vertex& _rhs) const {                             //加法
			return Vertex(x() + _rhs.x(), y() + _rhs.y(), z() + _rhs.z());
		}
		inline Vertex operator-(const Vertex& _rhs) const {                             //减法
			return Vertex(x() - _rhs.x(), y() - _rhs.y(), z() - _rhs.z());
		}
		inline Vertex operator*(double k) const {                                        //数乘
			return Vertex(x() * k, y() * k, z() * k);
		}
		inline Vertex operator/(double k) const {                                        //数除
			return Vertex(x() / k, y() / k, z() / k);
		}
		inline double operator*(const Vertex& _rhs) const {                              //点乘
			return double(x() * _rhs.x() + y() * _rhs.y() + z() * _rhs.z());
		}
		inline Vertex operator%(const Vertex& _rhs) const {                             //叉乘
			return Vertex(y() * _rhs.z() - z() * _rhs.y(),
				z() * _rhs.x() - x() * _rhs.z(),
				x() * _rhs.y() - y() * _rhs.x());
		}
		inline double norm() const {                                                      //模长
			return sqrt(x() * x() + y() * y() + z() * z());
		}
		inline Vertex normlize() const {                                                 //单位化
			auto m = this->norm();
			return Vertex(x() / m, y() / m, z() / m);
		}

		inline double dist(const Vertex& _rhs) const {                                   //计算俩个向量的距离
			return (*this - _rhs).norm();
		}


		/*===========比较运算=================*/
		inline bool operator==(const Vertex& _rhs) const {
			return (x() == _rhs.x() && y() == _rhs.y() && z() == _rhs.z());
		}
		/*===========赋值运算=================*/
		//返回值为引用，支持连等
		Vertex& operator=(const Vertex& _v) {
			Vec3 = _v.Vec3; return *this;
		}

		void setPosition(double x, double y, double z);

	private:
		std::vector<double> Vec3;
	};
}

/*======================特化V3f的哈希映射============================*/
//冲突时会调用相等函数
namespace std
{
	template<> struct hash<MeshKernel::Vertex>
	{
		size_t operator()(const MeshKernel::Vertex& v)const
		{
			return hash<double>()(v.x()) ^ hash<double>()(v.y()) ^ hash<double>()(v.z());
		}
	};
}