#pragma once
#include"Face.h"

namespace MeshKernel { 
	class Cell {
	public:
		/*=========================���캯��=============================*/
		Cell() :n_(0) {}
		Cell(const std::vector<VertexHandle>& _vertices, const std::vector<EdgeHandle>& _edges,
			 const std::vector<FaceHandle>& _faces){
			assert(_faces.size() + _vertices.size() - _edges.size() == 2);     //�Ƿ�����ŷ����ʽ
			assert(_faces.size() >= 4 && _vertices.size() >= 4);
			n_ = _faces.size();
			vertices_.assign(_vertices.begin(), _vertices.end());
			edges_.assign(_edges.begin(), _edges.end());
			faces_.assign(_faces.begin(), _faces.end());
		}
		Cell(const Cell& _c) {
			*this = _c;
		}
		/*=======================����������==========================*/
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
		/*=======================��дԪ��==========================*/
		const VertexHandle vh(int k) const { return vertices_[k]; }          // ��������ֻ֧�ֶ�
		const EdgeHandle eh(int k) const { return edges_[k]; }
		const FaceHandle fh(int k) const { return faces_[k]; }
		
		VertexHandle& vh(int k) { return vertices_[k]; }                     // �ǳ�������֧�ֶ�д
		EdgeHandle& eh(int k) { return edges_[k]; }
		FaceHandle& fh(int k) { return faces_[k]; }
		
		size_t size() const { return n_; }
		/*================�õ�����Ķ���handle����==========================*/
		std::vector<VertexHandle> getSortedVertexHandle() const {
			std::vector<VertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}
	private:
		// �ж�����������Ƿ��뵱ǰ��ȫ����ͬ
		bool isSameFaces(const Cell& _c) const {
			std::unordered_set<FaceHandle> fhset;
			for (const auto& f1 : _c.faces_) fhset.insert(f1);
			for (const auto& f2 : _c.faces_) fhset.insert(f2);
			return fhset.size() == n_;
		}
	private:
		// ����������е��handle
		std::vector<VertexHandle> vertices_;
		// ����������бߵ�handle
		std::vector<EdgeHandle> edges_;
		// ��������������handle
		std::vector<FaceHandle> faces_;
		// ������������
		size_t n_;                     
	};

}

/*======================�ػ�Cell�Ĺ�ϣӳ��============================*/
// To do: �޸�Hashӳ��
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