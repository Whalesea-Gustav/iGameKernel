#pragma once
#include"Edge.h"
#include<algorithm>
namespace MeshKernel{
	// �����е���Ƭ
	class Face {
	public: 
		/*=========================���캯��=============================*/
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
		/*=======================����������==========================*/
		inline Face& operator=(const Face& _f) {
			n_ = _f.n_;
			edges_.assign(_f.edges_.begin(), _f.edges_.end());
			vertices_.assign(_f.vertices_.begin(), _f.vertices_.end());
			return *this;
		}
		bool operator==(const Face& _f) const {
			return (n_== _f.n_) && isSameEdges(_f);
		}

		/*=======================��дԪ��==========================*/
		const VertexHandle vh(int k) const { return vertices_[k]; }            // ��������ֻ֧�ֶ�
		const EdgeHandle eh(int k) const { return edges_[k]; }
		VertexHandle& vh(int k) { return vertices_[k]; }		               // �ǳ�������֧�ֶ�д
		EdgeHandle& eh(int k) { return edges_[k]; }
		
		const size_t size() const { return n_; }

		/*================�õ�����Ķ���handle����==========================*/
		std::vector<VertexHandle> getSortedVertexHandle() const {
			std::vector<VertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}
	private:
		// �ж�������ı��Ƿ��뵱ǰ��ȫ����ͬ
		bool isSameEdges(const Face& _f) const {
			std::unordered_set<EdgeHandle> ehset;
			for (const auto& e1 : _f.edges_) ehset.insert(e1);
			for (const auto& e2 : edges_) ehset.insert(e2);
			return ehset.size() == n_;
		}
	private:
		// ����������бߵ�handle
		std::vector<EdgeHandle> edges_;
		// ����������е��handle
		std::vector<VertexHandle> vertices_;
		// �������ı���
		size_t n_;
	};	
}

/*======================�ػ�Face�Ĺ�ϣӳ��============================*/
// To do: �޸�Hashӳ��
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