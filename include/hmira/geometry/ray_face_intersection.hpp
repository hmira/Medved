#ifndef __LINE_FACE_INTERSECTION_H__
#define __LINE_FACE_INTERSECTION_H__

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <math.h>

template< typename Vector3T >
bool line_face_intersection(
	const Vector3T& rayStart,
	Vector3T& rayDir,
	const Vector3T& a,
	const Vector3T& b,
	const Vector3T& c,
	Vector3T& intersection,
	float& distance)
{
	
	distance = 0;
	intersection = Vector3T(0,0,0);   

	rayDir.normalize();

	Vector3T triNorm = cross(b-a, c-a).normalize();
	float vn = dot(rayDir, triNorm);
	
	if (vn < 0.0000001f && vn > -0.0000001f)
		return false;

	Vector3T aa = rayStart - a;
	float xpn = dot(aa, triNorm);
	distance = -xpn / vn;

	if (distance < 0) return false; // behind ray origin. fail

	auto hitPos = (rayDir * distance) + rayStart;

	Vector3T hit00 = hitPos - a;
	Vector3T hit01 = b - a;
	Vector3T cross0 = cross(hit00, hit01);
	if (dot(cross0, triNorm) > 0.0000001f) return false; // out of bounds

	Vector3T hit10 = hitPos - b;
	Vector3T hit11 = c - b;        
	Vector3T cross1 = cross(hit10, hit11);        
	if (dot(cross1, triNorm) > 0.0000001f) return false; // out of bounds

	Vector3T hit20 = hitPos - c;
	Vector3T hit21 = a - c;        
	Vector3T cross2 = cross(hit20, hit21);        
	if (dot(cross2, triNorm) > 0.0000001f) return false; // out of bounds
	
	intersection = hitPos;
	return true;
}

#endif //__LINE_FACE_INTERSECTION_H__