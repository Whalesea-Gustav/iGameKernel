#pragma once
#include<iostream>

namespace MeshKernel { 
	class Handle {
	public:
		// ��ʽ���죬������ "Handle h = 1;" ���﷨������ʽת��
		explicit Handle(int _idx) : idx_(_idx) {};
		// ��ֵ����
		Handle& operator=(const Handle& _h) {
			idx_ = _h.idx_;
			return *this;
		}

		// ����handle������int�����޸�
		Handle& operator=(int _idx) {
			idx_ = _idx;
			return *this;
		}

		//�ж�handle�Ƿ����
		inline bool is_valid() const { return idx_ != -1; }
		/*===========================================handle�ıȽϲ���===========================================*/
		inline bool operator<(const Handle& _h) const { return (this->idx_ < _h.idx_); }

		inline bool operator<(int _idx) const { return idx_ < _idx; }

		inline bool operator>(const Handle& _h) const { return (this->idx_ > _h.idx_); }

		inline bool operator>(int _idx) const { return idx_ > _idx; }

		inline bool operator==(const Handle& _h) const { return _h.idx_ == this->idx_; }

		inline bool operator!=(const Handle& _h) const { return _h.idx_ != this->idx_; }

		/*===========================================�޸�������===========================================*/

		inline const int& idx() const { return idx_; }      //ȡԪ��

		void idx(const int& _idx) { idx_ = _idx; }      //�޸�Ԫ��

		inline operator int() const { return idx_; } //��ʽת��Ϊint

		void reset() { idx_ = -1; }  //��ʼ��

	private:
		int idx_;
	};

	//��ʽ�����ܹ������������͵�handle��ֵ����ǰ���͵�handle
	class VertexHandle : public Handle { public: explicit VertexHandle(int _idx = -1) : Handle(_idx) {} };
	class EdgeHandle : public Handle { public: explicit EdgeHandle(int _idx = -1) : Handle(_idx) {} };
	class FaceHandle : public Handle { public: explicit FaceHandle(int _idx = -1) : Handle(_idx) {} };
	class CellHandle : public Handle { public: explicit CellHandle(int _idx = -1) : Handle(_idx) {} };
}

/*======================�ػ�����Handle�Ĺ�ϣӳ��============================*/
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
// �������������û����ѱ����ǵĲ�����
// ���� "Handle& operator=(int _idx)"
//void test() {
//	VertexHandle vh(1);
//	vh = 2;                      ���󣬲�������ʽת��
//	vh.Handle::operator= (2);    ��ȷ�����û��������
//}