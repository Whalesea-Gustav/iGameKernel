#include "Vertex.h"

namespace  MeshKernel {
	void Vertex::setPosition(double x, double y, double z)
	{	
		this->Vec3[0] = x;
		this->Vec3[1] = y;
		this->Vec3[2] = z;
	}
}
