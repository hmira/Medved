#ifndef __LINE_ON_PLANE_PROJECTION_HPP__
#define __LINE_ON_PLANE_PROJECTION_HPP__

#include <iostream>
#include <limits>
#include <vector>
#include <math.h>
#include <boost/concept_check.hpp>
#include <boost/config/no_tr1/complex.hpp>

namespace hmira
{

namespace geometry
{
	template <typename TVector>
	void two_lines_intersection(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2);
	
	template <typename TVector>
	void two_lines_intersection_x(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2);
	
	template <typename TVector>
	void two_lines_intersection_y(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2);
	
	template <typename TVector>
	void two_lines_intersection_z(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2);
	
	template <typename TFace>
	void line_on_face_projection(TFace face)
	{
		
	}
	
	template <typename TVector>
	bool
	line_on_plane_projection(
		const TVector& a,
		const TVector& b,
		const TVector& c,
		const TVector& face_normal,
		const TVector& line_origin,
		const TVector& line_dir,
		TVector& proj_point1,
		TVector& proj_point2
				)
	{
		std::vector<TVector> result;
		
		auto dir = line_dir;
		dir.normalize();
		
		auto normal = face_normal;
		normal.normalize();
		
		if (dot(dir, normal) == 1 || dot(dir, normal) == -1)
		{
			proj_point1 = line_origin;
			proj_point2 = line_origin;
			return true;
		}
		
		auto aligned = cross(dir, normal);
		
		auto proj = cross(normal, aligned);
		
		float t1_ab = std::numeric_limits<float>::infinity();
		float t2_ab = std::numeric_limits<float>::infinity();
		two_lines_intersection(a, (b-a), line_origin, proj, t1_ab, t2_ab);
		if ((t1_ab >= 0 && t1_ab <= 1))
			result.push_back(a + t1_ab*(b-a));
		
		float t1_ac = std::numeric_limits<float>::infinity();
		float t2_ac = std::numeric_limits<float>::infinity();
		two_lines_intersection(a, (c-a), line_origin, proj, t1_ac, t2_ac);
		if ((t1_ac >= 0 && t1_ac <= 1))
			result.push_back(a + t1_ac*(c-a));
		
		float t1_bc = std::numeric_limits<float>::infinity();
		float t2_bc = std::numeric_limits<float>::infinity();
		two_lines_intersection(c, (b-c), line_origin, proj, t1_bc, t2_bc);
		if ((t1_bc >= 0 && t1_bc <= 1) && (result.size() == 1))
			result.push_back(c + t1_bc*(b-c));
		
		if (result.size() == 2)
		{
			proj_point1 = result[0];
			proj_point2 = result[1];
			return true;
		}
		else
		{
			return false;
		}
	}
	
	template <typename TVector>
	void two_lines_intersection(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2)
	{
		two_lines_intersection_x(p1, v1, p2, v2, t1, t2);
		if (!isnan(t1)) return;
		two_lines_intersection_y(p1, v1, p2, v2, t1, t2);
		if (!isnan(t1)) return;
		two_lines_intersection_z(p1, v1, p2, v2, t1, t2);
		if (!isnan(t1)) return;
	}

	/**
	 * 
	 * @return point in 3D/2D
	 */
	template <typename TVector>
	void two_lines_intersection_x(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2)
	{
		t1 = (v2[1]*p1[2] - v2[1]*p2[2] - v2[2]*p1[1] + v2[2]*p2[1]) / (v2[2]*v1[1] - v2[1]*v1[2]);
		t2 = (p1[1] + t1*v1[1] - p2[1]) / v2[1];
		
// 		return p1 + t1 * v1;
	}
	
	
	/**
	 * 
	 * @return point in 3D/2D
	 */
	template <typename TVector>
	void two_lines_intersection_y(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2)
	{
		t1 = (v2[0]*p1[2] - v2[0]*p2[2] - v2[2]*p1[0] + v2[2]*p2[0]) / (v2[2]*v1[0] - v2[0]*v1[2]);
		t2 = (p1[0] + t1*v1[0] - p2[0]) / v2[0];
		
// 		return p1 + t1 * v1;
	}
	
	
	/**
	 * 
	 * @return point in 3D/2D
	 */
	template <typename TVector>
	void two_lines_intersection_z(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2,
		float& t1,
		float& t2)
	{
		t1 = (v2[0]*p1[1] - v2[0]*p2[1] - v2[1]*p1[0] + v2[1]*p2[0]) / (v2[1]*v1[0] - v2[0]*v1[1]);
		t2 = (p1[0] + t1*v1[0] - p2[0]) / v2[0];
		
// 		return p1 + t1 * v1;
	}
	
	
	/**
	 * 
	 * @return point in 3D/2D
	 */
	template <typename TVector>
	TVector two_lines_intersection(
		const TVector& p1,
		const TVector& v1,
		const TVector& p2,
		const TVector& v2)
	{
		auto s = (v2[0]*p1[1] - v2[0]*p2[1] - v2[1]*p1[0] + v2[1]*p2[0]) / (v2[1]*v1[0] - v2[0]*v1[1]);
		auto t = (p1[0] + s*v1[0] - p2[0]) / v2[0];
		
		return p1 + s * v1;
	}
	
} //geometry

} //hmira

#endif //__LINE_ON_PLANE_PROJECTION_HPP__