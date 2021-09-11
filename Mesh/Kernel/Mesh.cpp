#include"Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
// Mesh 定义
namespace MeshKernel { 
	

	FaceHandle Mesh::AddFace(const std::vector<VertexHandle>& _vhs) {
		std::vector<EdgeHandle> ehs(_vhs.size());
		// 得到该面边的handle
		for (int i = 0; i < _vhs.size(); ++i) {
			if (i == 0) {
				ehs[i] = AddEdge(_vhs[_vhs.size() - 1], _vhs[i]);
			}
			else {
				ehs[i] = AddEdge(_vhs[i], _vhs[i - 1]);
			}
		}
		Face f(_vhs, ehs);
		// 如果该面已经存在，则返回面的handle
		if (Face2Fh_.count(f)) {
			return Face2Fh_[f];
		}
		// 否则
		else {
			FaceHandle fh = GenFaceHandle();						// 生成一个新的handle
			faces_[fh] = f;											// 建立handle与该面之间的映射
			Face2Fh_[f] = fh;										// 建立该面与handle之间的映射
			AddFace2Neighbor(fh);									// 将该面添加至面所包含的点和边的相邻面中
			return Face2Fh_[f];										// 返回该面的handle
		}
	}
	EdgeHandle Mesh::AddEdge(const VertexHandle& vh1, const VertexHandle& vh2) {
		Edge e(vh1, vh2);
		// 如果该边已经存在，则返回该边的handle
		if (Edge2Eh_.count(e)) return Edge2Eh_[e];
		// 否则
		else {
			EdgeHandle eh = GenEdgeHandle();						// 生成一个新的handle
			edges_[eh] = e;											// 建立handle与该边之间的映射
			Edge2Eh_[e] = eh;										// 建立该边与handle之间的映射
			AddEdge2Neighbor(eh);									// 将该边记录为其两个顶点的邻接边
			return Edge2Eh_[e];										// 返回该边的handle
		}
	}
	VertexHandle Mesh::AddVertex(const Vertex& _v) {
		
		if (Vertex2Vh_.count(_v)) return Vertex2Vh_[_v];
		else {
			VertexHandle vh = GenVertexHandle();
			vertices_[vh] = _v;
			Vertex2Vh_[_v] = vh;
			return Vertex2Vh_[_v];
		}
	
	}
	VertexHandle Mesh::DeleteVertex(const VertexHandle& _vh) {
		if (!vertices_.count(_vh)) return VertexHandle(-1);        // 如果该点不存在，则返回-1
		else {
			// 删除相邻边元素（删除边时自然删除了相邻面）
			auto ve = NeighborEhOfVertex_[_vh];
			for (EdgeHandle eh : ve) 
			{
				Edge e = this->edges(eh);
				//std::cout << "delete edge with vertices :  " << e.vh1() << " | " << e.vh2() << std::endl;

				DeleteEdge(eh);
				
			}
			// 删除顶点元素和以其为根据的所有邻接关系
			vertices_.erase(_vh);
			NeighborEhOfVertex_.erase(_vh);
			NeighborFhOfVertex_.erase(_vh);
			return _vh;
		}

	}
	EdgeHandle Mesh::DeleteEdge(const EdgeHandle& _eh) {
		if (!edges_.count(_eh)) return EdgeHandle(-1);             // 如果该边不存在，则返回-1
		else {
			// 删除邻接关系
			Edge e(edges_[_eh]);
			for (int i = 0; i < 2; ++i) {
				VertexHandle ev = e.vh(i);
				NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
			}
			// 删除相邻面元素
			auto ef = NeighborFhOfEdge_[_eh];
			for (FaceHandle fh : ef ) {
				DeleteFace(fh);
			}
			// 删除边元素和以其为根据的所有邻接关系
			edges_.erase(_eh);
			NeighborFhOfEdge_.erase(_eh);
			return _eh;
		}
	}
	FaceHandle Mesh::DeleteFace(const FaceHandle& _fh) {
		if (!faces_.count(_fh)) return FaceHandle(-1);             // 如果该面不存在，则返回-1
		else {                                                     // 如果该面存在，则返回删除的这个面的handle
			// 删除邻接关系
			Face f(faces_[_fh]);
			for (int i = 0; i < f.size(); ++i) {
				VertexHandle fv = f.vh(i);
				EdgeHandle fe = f.eh(i);
				NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
				NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
			}
			// 删除面元素
			faces_.erase(_fh);
			return _fh;
		}
	}

	void Mesh::AddFace2Neighbor(const FaceHandle& _fh)
	{
		Face f = faces_[_fh];
		size_t n = f.size();
		for (int i = 0; i < n; ++i) {
			NeighborFhOfVertex_[f.vh(i)].insert(_fh);
		}
		for (int i = 0; i < n; ++i) {
			NeighborFhOfEdge_[f.eh(i)].insert(_fh);
		}
	}
	void Mesh::AddEdge2Neighbor(const EdgeHandle& _eh)
	{
		Edge e = edges_[_eh];
		NeighborEhOfVertex_[e.vh1()].insert(_eh);
		NeighborEhOfVertex_[e.vh2()].insert(_eh);

	}
	void Mesh::DeleteFace2Neighbor(const FaceHandle& _fh){
		Face f = faces_[_fh];
		size_t n = f.size();
		for (int i = 0; i < n; ++i) {
			NeighborFhOfVertex_[f.vh(i)].erase(_fh);
		}
		for (int i = 0; i < n; ++i) {
			NeighborFhOfEdge_[f.eh(i)].erase(_fh);
		}
	}
	void Mesh::DeleteEdge2Neighbor(const EdgeHandle& _eh) {
		Edge e = edges_[_eh];
		NeighborEhOfVertex_[e.vh1()].erase(_eh);
		NeighborEhOfVertex_[e.vh2()].erase(_eh);
	}

	Mesh & Mesh::operator=(const Mesh & _surfacemesh)
	{
		vertices_ = _surfacemesh.vertices_;
		edges_ = _surfacemesh.edges_;
		faces_ = _surfacemesh.faces_;
		Vertex2Vh_ = _surfacemesh.Vertex2Vh_;
		Edge2Eh_ = _surfacemesh.Edge2Eh_;
		Face2Fh_ = _surfacemesh.Face2Fh_;
		NeighborEhOfVertex_ = _surfacemesh.NeighborEhOfVertex_;
		NeighborFhOfVertex_ = _surfacemesh.NeighborFhOfVertex_;
		NeighborFhOfEdge_ = _surfacemesh.NeighborFhOfEdge_;
		return *this;

	}

	/*=========================读写元素===============================*/
	// 读取ID为i的顶点
	Vertex& Mesh::vertices(VertexHandle _vh) {
		assert(vertices_.count(_vh));
		return vertices_[_vh];
	}
	const Vertex Mesh::vertices(VertexHandle _vh) const {
		assert(vertices_.count(_vh));
		return vertices_.find(_vh)->second;                // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用   
	}
	// 读取ID为i的边
	Edge& Mesh::edges(EdgeHandle _eh) {
		assert(edges_.count(_eh));
		return edges_[_eh];
	}
	const Edge& Mesh::edges(EdgeHandle _eh) const {
		assert(edges_.count(_eh));
		return edges_.find(_eh)->second;
	}
	// 读取ID为i的面
	Face& Mesh::faces(FaceHandle _fh) {
		assert(faces_.count(_fh));
		return faces_[_fh];
	}
	const Face Mesh::faces(FaceHandle _fh) const {
		assert(faces_.count(_fh));
		return faces_.find(_fh)->second;
	}

	/*====================根据元素得到对应ID=========================*/
	const VertexHandle Mesh::vertexhanle(Vertex _vertex) const {
		if (Vertex2Vh_.find(_vertex) != Vertex2Vh_.end()) return Vertex2Vh_.find(_vertex)->second;
		else return VertexHandle(-1);
	}
	const EdgeHandle Mesh::edgehandle(Edge& _edge) const {
		if (Edge2Eh_.find(_edge) != Edge2Eh_.end()) return Edge2Eh_.find(_edge)->second;
		else return EdgeHandle(-1);
	}
	const FaceHandle Mesh::facehandle(Face& _face) const {
		if (Face2Fh_.find(_face) != Face2Fh_.end()) return Face2Fh_.find(_face)->second;
		else return FaceHandle(-1);
	}

	/*======================得到邻接关系============================*/
	// 顶点的邻接点
	// 先找邻接边，再找邻接点
	std::unordered_set<VertexHandle> Mesh::NeighborVh(VertexHandle _vh) {
		std::unordered_set<VertexHandle> neighborvh;
		auto neighboreh = NeighborEh(_vh);
		// 存在邻接边是前提
		if (neighboreh.size()) {
			for (EdgeHandle eh : neighboreh) {
				if (edges_[eh].vh1() != _vh) neighborvh.insert(edges_[eh].vh1());
				if (edges_[eh].vh2() != _vh) neighborvh.insert(edges_[eh].vh2());
			}
		}
		return neighborvh;
	}
	// 顶点的邻接边
	std::unordered_set<EdgeHandle> Mesh::NeighborEh(VertexHandle _vh) {
		if (NeighborEhOfVertex_.count(_vh)) return NeighborEhOfVertex_[_vh];
		else return std::unordered_set<EdgeHandle>();               // 返回一个空的集合
	}
	// 顶点的邻接面
	std::unordered_set<FaceHandle> Mesh::NeighborFh(VertexHandle _vh) {
		if (NeighborFhOfVertex_.count(_vh)) return NeighborFhOfVertex_[_vh];
		else return std::unordered_set<FaceHandle>();               // 返回一个空的集合
	}
	// 边的邻接边
	// 两个顶点的所有邻接边去除当前边
	std::unordered_set<EdgeHandle> Mesh::NeighborEh(EdgeHandle _eh) {
		assert(edges_.count(_eh));                     // 保证该边handle存在
		std::unordered_set<EdgeHandle> neighboreh;     // 保存输出的结果
		int k = 0;                                     // 遍历两个顶点
		while (k < 2) {                                
			VertexHandle vh = edges_[_eh].vh(k);
			auto vhneighboreh = NeighborEh(vh);        // 得到点的邻接边
			for (EdgeHandle eh : vhneighboreh) {
				if (eh != _eh) neighboreh.insert(eh);
			}
			++k;                  
		}
		
		return neighboreh;
	}
	// 边的邻接面
	std::unordered_set<FaceHandle> Mesh::NeighborFh(EdgeHandle _eh) {
		if (NeighborFhOfEdge_.count(_eh)) return NeighborFhOfEdge_[_eh];
		else return std::unordered_set<FaceHandle>();               // 返回一个空的集合
	}
	// 面的邻接面
	// 邻接面：有一条相同边
	std::unordered_set<FaceHandle> Mesh::NeighborFh(FaceHandle _fh) {
		assert(faces_.count(_fh));                     // 保证该边handle存在
		std::unordered_set<FaceHandle> neigborface;
		int k = 0;                                     // 遍历两个顶点
		size_t facesize = faces_[_fh].size();
		while (k < facesize) {
			EdgeHandle eh = faces_[_fh].eh(k);
			auto ehneighborfh = NeighborFh(eh);        // 得到点的邻接边
			for (FaceHandle fh : ehneighborfh) {
				if (fh != _fh) neigborface.insert(fh);
			}
			++k;
		}
		return neigborface;
	}
}



// SurfaceMesh 定义
namespace MeshKernel {
	void SurfaceMesh::InitMesh(const std::vector<Vertex>& _vertices,
		const std::vector<std::vector<VertexHandle>>& _elements) {
		for (auto v : _vertices) {
			auto vh = AddVertex(Vertex(v.x(), v.y(), v.z()));
		}
		for (auto f : _elements) {
			AddFace(f);
		}
	}
	SurfaceMesh& SurfaceMesh::operator=(const SurfaceMesh& _surfacemesh) {
		if (this != &_surfacemesh) {
			Mesh::operator=(_surfacemesh);
		}
		return *this;
	}

	bool SurfaceMesh::IsOnBoundary(VertexHandle vh)
	{
		std::unordered_set<EdgeHandle> neighEh = this->NeighborEh(vh);
		for (EdgeHandle eh : neighEh)
		{
			if (this->IsOnBoundary(eh))
			{
				return true;
			}
		}
		return false;
	}
	bool SurfaceMesh::IsOnBoundary(EdgeHandle eh)
	{
		return this->NeighborFh(eh).size() < 2;
	}

	void SurfaceMesh::collapse(EdgeHandle& eh)
	{
		Edge e = this->edges(eh);
		VertexHandle vh1 = e.vh1();
		VertexHandle vh2 = e.vh2();
		std::cout << vh1.idx() << " | " << vh2.idx() << std::endl;

		std::cout << "before adding, start testing for vh2" << std::endl;
		for (FaceHandle neighborFh : this->NeighborFh(vh2))
		{
			Face neighborFace = this->faces(neighborFh);
			std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		}

		for (FaceHandle neighborFh : this->NeighborFh(vh1))
		{

			Face face = this->faces(neighborFh);
			VertexHandle& fvh0 = face.vh(0);
			VertexHandle& fvh1 = face.vh(1);
			VertexHandle& fvh2 = face.vh(2);
			//std::cout << fvh0.idx() << " | " << fvh1.idx() << " | " << fvh2.idx() << std::endl;

			if (fvh0 == vh2 || fvh1 == vh2 || fvh2 == vh2)
			{
				//包含vh1, vh2的面需要被删除
				continue;
			}

			std::vector<VertexHandle> MyFace(3);

			if (fvh0 == vh1)
			{
				//wrong order influence rendering (keep counterclock direction)
				MyFace[0] = vh2;
				MyFace[1] = fvh1;
				MyFace[2] = fvh2;
				//std::cout << vh2.idx() << " " << fvh1.idx() << " " << fvh2.idx() << std::endl;
			}
			else if(fvh1 == vh1)
			{
				MyFace[0] = fvh0;
				MyFace[1] = vh2;
				MyFace[2] = fvh2;
				//std::cout << vh2.idx() << " " << fvh0.idx() << " " << fvh2.idx() << std::endl;
			}
			else if (fvh2 == vh1)
			{
				MyFace[0] = fvh0;
				MyFace[1] = fvh1;
				MyFace[2] = vh2;
				//std::cout << vh2.idx() << " " << fvh0.idx() << " " << fvh1.idx() << std::endl;
			}

			std::cout << MyFace[0].idx() << " | " << MyFace[1].idx() << " | " << MyFace[2].idx() << std::endl;
			
			AddFace(MyFace);
		}

		//test for collapse edge 90

		std::cout << "start testing for vh2" << std::endl;
		for (FaceHandle neighborFh : this->NeighborFh(vh2))
		{
			Face neighborFace = this->faces(neighborFh);
			std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		}

		DeleteVertex(vh1);

		//std::cout << vh2.idx() << std::endl;
		std::cout << "after delete vh1, start testing for vh2" << std::endl;
		for (FaceHandle neighborFh : this->NeighborFh(vh2))
		{
			Face neighborFace = this->faces(neighborFh);
			std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		}
	}
	void SurfaceMesh::collapse2(EdgeHandle& eh)
	{
		Edge& e = this->edges(eh);
		VertexHandle& vh1 = e.vh1();
		VertexHandle& vh2 = e.vh2();
		Vertex v2 = this->vertices(vh2);
		Vertex v1 = this->vertices(vh1);
		std::cout << vh1.idx() << " | " << vh2.idx() << std::endl;

		Vertex newVertex(0.5 *(v1.x() +  v2.x()), 0.5 * (v1.y() + v2.y()), 0.5 * (v1.z() + v2.z()));
		VertexHandle newVh = this->AddVertex(newVertex);
		std::cout << "new vh idx : " << newVh.idx() << std::endl;

		std::vector<EdgeHandle> one_ring_edge;

		//查找vh1的邻面
		for (FaceHandle neighborFh : this->NeighborFh(vh1))
		{
			
			Face& face = this->faces(neighborFh);
			std::cout << "face components : " << face.vh(0) << " | " << face.vh(1) << " | " << face.vh(2) << std::endl;
			for (int i = 0; i < 3; i++)
			{
				EdgeHandle& eh = face.eh(i);
				Edge& edge = this->edges(eh);
				if (edge.vh1() == vh1 || edge.vh1() == vh2 || edge.vh2() == vh1 || edge.vh2() == vh2)
				{
					continue;
				}
				else
				{
					one_ring_edge.push_back(eh);
				}
			}
		}
		//查找vh2的邻面
		for (FaceHandle neighborFh : this->NeighborFh(vh2))
		{
			Face& face = this->faces(neighborFh);
			std::cout << "face components : " << face.vh(0) << " | " << face.vh(1) << " | " << face.vh(2) << std::endl;
			
			for (int i = 0; i < 3; i++)
			{
				EdgeHandle& eh = face.eh(i);
				Edge& edge = this->edges(eh);

				if (edge.vh1() == vh1 || edge.vh1() == vh2 || edge.vh2() == vh1 || edge.vh2() == vh2)
				{
					continue;
				}
				else
				{
					one_ring_edge.push_back(eh);
				}
			}
		}

		for (EdgeHandle one_ring_eh : one_ring_edge)
		{
			std::vector<VertexHandle> newFace(3);
			Edge one_ring_e = this->edges(one_ring_eh);
			VertexHandle one_ring_vh1 = one_ring_e.vh1();
			VertexHandle one_ring_vh2 = one_ring_e.vh2();
			newFace[0] = newVh;
			newFace[1] = one_ring_vh1;
			newFace[2] = one_ring_vh2;
			std::cout << newFace[0].idx() << " | " << newFace[1].idx() << " | " << newFace[2].idx() << std::endl;
			AddFace(newFace);
		}
		DeleteVertex(vh1);
		DeleteVertex(vh2);
		
	}

	void SurfaceMesh::collapseTo(EdgeHandle& eh, double x, double y, double z)
	{
		Edge e = this->edges(eh);
		VertexHandle vh1 = e.vh1();
		VertexHandle vh2 = e.vh2();
		Vertex& v2 = this->vertices(vh2);

		v2.x() = x;
		v2.y() = y; 
		v2.z() = z;

		//std::cout << vh1.idx() << " | " << vh2.idx() << std::endl;

		//std::cout << "before adding, start testing for vh2" << std::endl;
		//for (FaceHandle neighborFh : this->NeighborFh(vh2))
		//{
		//	Face neighborFace = this->faces(neighborFh);
		//	std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		//}

		for (FaceHandle neighborFh : this->NeighborFh(vh1))
		{

			Face face = this->faces(neighborFh);
			VertexHandle& fvh0 = face.vh(0);
			VertexHandle& fvh1 = face.vh(1);
			VertexHandle& fvh2 = face.vh(2);
			//std::cout << fvh0.idx() << " | " << fvh1.idx() << " | " << fvh2.idx() << std::endl;

			if (fvh0 == vh2 || fvh1 == vh2 || fvh2 == vh2)
			{
				//包含vh1, vh2的面需要被删除
				continue;
			}

			std::vector<VertexHandle> MyFace(3);

			if (fvh0 == vh1)
			{
				//wrong order influence rendering (keep counterclock direction)
				MyFace[0] = vh2;
				MyFace[1] = fvh1;
				MyFace[2] = fvh2;
				//std::cout << vh2.idx() << " " << fvh1.idx() << " " << fvh2.idx() << std::endl;
			}
			else if (fvh1 == vh1)
			{
				MyFace[0] = fvh0;
				MyFace[1] = vh2;
				MyFace[2] = fvh2;
				//std::cout << vh2.idx() << " " << fvh0.idx() << " " << fvh2.idx() << std::endl;
			}
			else if (fvh2 == vh1)
			{
				MyFace[0] = fvh0;
				MyFace[1] = fvh1;
				MyFace[2] = vh2;
				//std::cout << vh2.idx() << " " << fvh0.idx() << " " << fvh1.idx() << std::endl;
			}

			//std::cout << MyFace[0].idx() << " | " << MyFace[1].idx() << " | " << MyFace[2].idx() << std::endl;

			AddFace(MyFace);
		}

		//test for collapse edge 90

		//std::cout << "start testing for vh2" << std::endl;
		//for (FaceHandle neighborFh : this->NeighborFh(vh2))
		//{
		//	Face neighborFace = this->faces(neighborFh);
		//	std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		//}

		DeleteVertex(vh1);

		////std::cout << vh2.idx() << std::endl;
		//std::cout << "after delete vh1, start testing for vh2" << std::endl;
		//for (FaceHandle neighborFh : this->NeighborFh(vh2))
		//{
		//	Face neighborFace = this->faces(neighborFh);
		//	std::cout << "Face " << neighborFh.idx() << " : " << neighborFace.vh(0).idx() << " | " << neighborFace.vh(1).idx() << " | " << neighborFace.vh(2).idx() << std::endl;
		//}
	}

	void SurfaceMesh::test_iteration()
	{
		VertexHandle vh(695);
		for (auto fh : this->NeighborFh(vh))
		{
			std::cout << fh.idx() << std::endl;
		}

		//启示: Face.vh(i)保证逆时针有序，Face.eh(i)只提供连接关系，并不能保证有序
		for (auto fh : this->NeighborFh(vh))
		{
			Face f = this->faces(fh);
			std::cout << "iter face id : " << fh.idx() << std::endl;
			for (int j = 0; j < 3; j++)
			{
				EdgeHandle& eh_f = f.eh(j);
				Edge& e_f = this->edges(eh_f);
				std::cout << "edges of face : " << e_f.vh1().idx() << " and " << e_f.vh2().idx() << std::endl;
			}
		}
	}

	//O(degree);
	EdgeHandle& SurfaceMesh::findEdge(VertexHandle& from_v, VertexHandle& to_v)
	{
		for (EdgeHandle neighborEh : this->NeighborEh(from_v))
		{
			Edge& e = this->edges(neighborEh);

			if (e.vh1() == to_v || e.vh2() == to_v)
			{
				return neighborEh;
			}
		}
		return EdgeHandle(-1);
	}

	//EdgeHandle(90).vh1()
	//use vertexhandle(695) as test exaple
	std::vector<EdgeHandle> SurfaceMesh::one_ring_edges(VertexHandle& vh)
	{
		std::unordered_set<FaceHandle> neighborFh = this->NeighborFh(vh);
		
		std::vector<EdgeHandle> results_in_order;
		
		if (neighborFh.size() == 0)	return results_in_order;

		//方法1：找到所有的对边，首尾相连
		//方法2:	 依顺序找到所有的对边

		//bool firstIter = true;
		//VertexHandle preVh(-1);

		//for (auto fh : this->NeighborFh(vh))
		//{
		//	std::cout << fh.idx() << std::endl;
		//}

		/*for (auto iter = neighborFh.begin(); iter != neighborFh.end(); ++iter)
		{
			if (firstIter)
			{
				firstIter = false;
				FaceHandle fh = *iter;
				Face f = this->faces(fh);
				for (int i = 0; i < 3; i++)
				{
					EdgeHandle eh_f = f.eh(i);
					Edge e_f = this->edges(eh_f);
					if (e_f.vh1() == vh || e_f.vh2() == vh)
					{
						continue;
					} 
					else
					{
						std::cout << "Face idx : " << fh.idx() << " | " << "vertex idx : " << e_f.vh1().idx() << std::endl;
						preVh = e_f.vh2();
						results_in_order.push_back(eh_f);					
					}
				}
				break;
			}
		}
		
		for (auto fh : this->NeighborFh(vh))
		{
			std::cout << fh.idx() << std::endl;
		}*/

		//for (int i = 0; i < 1; i++)
		//{
		//	for (auto iter = neighborFh.begin(); iter != neighborFh.end(); ++iter)
		//	{
		//		FaceHandle fh = *iter;
		//		Face f = this->faces(fh);
		//		for (int j = 0; j < 3; j++)
		//		{
		//			EdgeHandle eh_f = f.eh(j);
		//			Edge e_f = this->edges(eh_f);
		//			if (e_f.vh1() == preVh)
		//			{
		//				std::cout << "Face idx : " << fh.idx() << " | " << "vertex idx : " << e_f.vh1().idx() << std::endl;
		//				preVh = e_f.vh2();
		//				results_in_order.push_back(eh_f);		
		//			}
		//		}
		//	}
		//}
	}

	std::vector<EdgeHandle> SurfaceMesh::get_one_ring_edges(VertexHandle& vh)
	{
		std::vector<EdgeHandle> results;
		if (this->NeighborFh(vh).size() == 0)
		{
			return results;
		}
		bool firstIter = true;
		VertexHandle preVh;

		//方法1: 模仿半边结构，依顺序找到所有的one_ring_edges
		for (auto fh : this->NeighborFh(vh))
		{
			//std::cout << "iter face id : " << fh.idx() << std::endl;
			Face f = this->faces(fh);
			if (firstIter)
			{
				firstIter = false;
				
				for (int i = 0; i < 3; i++)
				{
					if (f.vh(i) == vh)
					{
						//find edge that does`t embrace vh
						preVh = f.vh((i + 2) % 3);
						VertexHandle from_vh = f.vh((i + 1) % 3);
						//std::cout << "edges of face : " <<from_vh.idx() << " and " << preVh.idx() << std::endl;

						Edge e = Edge(from_vh, preVh);
						results.push_back(this->edgehandle(e));
						break;
					}
				}
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					//find edge starting from preVh
					if (f.vh(i) == preVh)
					{
						preVh = f.vh((i + 1) % 3);
						//std::cout << "edges of face : " << f.vh(i).idx() << " and " << preVh.idx() << std::endl;
						Edge e = Edge(f.vh(i), preVh);
						results.push_back(this->edgehandle(e));
						break;
					}
				}
			}
		}
		
		int left_edge_count = this->NeighborFh(vh).size() - results.size();

		//std::cout << "iteration for rest edges." << std::endl;
		
		while (left_edge_count > 0)
		{
			for (auto fh : this->NeighborFh(vh))
			{
				//如果一次fh遍历就找完了所有的边，直接结束遍历
				if (left_edge_count == 0) break;

				Face f = this->faces(fh);
				//std::cout << "iter face id : " << fh.idx() << std::endl;
				for (int j = 0; j < 3; j++)
				{
					//find edge starting from preVh
					if (f.vh(j) == preVh)
					{
						preVh = f.vh((j + 1) % 3);
						//找到一条边
						left_edge_count--;
						//std::cout << "edges of face : " << f.vh(j).idx() << " and " << preVh.idx() << std::endl;
						Edge e = Edge(f.vh(j), preVh);
						results.push_back(this->edgehandle(e));
						break;
					}
				}
			}
		}
	
		return results;
	}

	std::vector<EdgeHandle> SurfaceMesh::get_one_ring_edges2(VertexHandle& vh)
	{
		std::vector<EdgeHandle> results;
		if (this->NeighborFh(vh).size() == 0)
		{
			return results;
		}
		//方法2: 无序地找到所有的 one_ring_edges，实现首尾相连算法

		//O(n)
		VertexHandle preVh;
		for (auto fh : this->NeighborFh(vh))
		{
			Face f = this->faces(fh);
			for (int i = 0; i < 3; i++)
			{
				if (f.vh(i) == vh)
				{				
					VertexHandle from_vh = f.vh((i + 1) % 3);
					if (preVh.idx() == -1)
					{
						preVh = from_vh;
					}
					VertexHandle to_vh = f.vh((i + 2) % 3);
					Edge e = Edge(from_vh, to_vh);
					results.push_back(this->edgehandle(e));
					break;
				}
			}
		}
		//to do 
		//1.找到一个起始点
		//2.首尾连接成环
		std::vector<bool> isAdded(results.size(),false);
		std::vector<EdgeHandle> ordered;

		while (ordered.size() < results.size())
		{
			for (int i = 0; i < results.size(); i++)
			{
				if (isAdded[i] == true) continue; //跳过已经添加了的边
				
				EdgeHandle& eh = results[i];
				Edge& e = this->edges(eh);
				
				if (e.vh1() == preVh)
				{
					ordered.push_back(eh);
					preVh = e.vh2();
					isAdded[i] = true;
				}
				else if (e.vh2() == preVh)
				{
					ordered.push_back(eh);
					preVh = e.vh1();
					isAdded[i] = true;
				}

			}
		}

		return ordered;
	}
}

//to do 
// 1. one-ring的对边遍历

