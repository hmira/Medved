#ifndef __POINT_IN_POLYHEDRON_HPP__
#define __POINT_IN_POLYHEDRON_HPP__

#include <iostream>
#include <unordered_map>
#include <tuple>

#include <hmira/range/faces.hpp>
#include <hmira/geometry/line_on_plane_projection.hpp>

#include <hmira/preprocessor/debug_mode.hpp>

namespace hmira
{

namespace geometry
{
	/**
	 * \brief possible cases of point and polyhedron in space
	 */
	enum point_classification
	{
		outside = 0,
		inside = 1,
		boundary = 2
	};
	
	/**
	 * \brief method that determines whether the point is either inside or outside the polygonal mesh
	 * 
	 * the method determines the point classification by counting the number of intersections
	 * led from the point in given direction. The direction is set by default to [1,1,1]
	 * 
	 * if the 
	 * 
	 * \return	the point classification - inside/outside/boundary
	 * 
	 * \param m	polyhedron represented by mesh
	 * 
	 */
	template <typename TPoint, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	point_classification point_in_polyhedron(
		TMesh& m,
		TPoint tested_point,
		TPoint tested_direction = TPoint(1,1,1))// = TMesh_Traits::Point(1, 1, 1))
	{
		std::unordered_map<float, int> distances; /// set of points that intersect the line led from the tested_point

		range::for_each_face(m,
			[&]
			(typename TMesh_Traits::face_descriptor face)
			{
				auto fv_pair = TMesh_Traits::get_surrounding_vertices(m, face);
				typename TMesh_Traits::vertex_descriptor v[3];
				
				float dist = 0.f;
				OpenMesh::Vec3f intersect(0,0,0);
				
				auto fv = fv_pair.first;
				for(int i = 0; i < 3; i++ )
				{
					v[i] = *fv;
					++fv;
				}
				
				auto a = TMesh_Traits::get_coordinates(m, v[0]);
				auto b = TMesh_Traits::get_coordinates(m, v[1]);
				auto c = TMesh_Traits::get_coordinates(m, v[2]);
			
				auto result = line_face_intersection(
				tested_point,
				tested_direction,
				a,
				b,
				c,
				intersect,
				dist);

				if (result)
				{
					
					auto dot_prod = dot(TMesh_Traits::get_face_normal(m, face), tested_direction) > 0;
					int entry_flag = 0; 
					if (dot_prod)
					{
						entry_flag = 1;
					}
					else
					{
						entry_flag = 2;
					}
					
					distances[dist] |= entry_flag;
					
				}
			});
		
		auto boundary = false;
		int number_of_intersects = 0;
		
		for (auto a:distances)
		{
			if ( a.first <= 0.0001f &&  a.first >= -0.0001f )
				boundary = true;

			if ( a.second == 1 || a.second == 2 )
				number_of_intersects++;
		}
		
		
		if (boundary)
		{
			
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the point is evaluated as boundary" )
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the position of the point refers to" )
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the surface of the object")
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
			return point_classification::boundary;
		}

		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the number of intersections with" )
		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : object and the line led from the" )
		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : given point: ", tested_point ", in direction: ", tested_direction )
		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : is ", number_of_intersects ," in total" )
		H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
		
		if ( number_of_intersects & 1 )
		{
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the point is inside of the object" )
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
			return point_classification::inside;
		}
		else
		{
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : the point is outside of the object" )
			H_DEBUG_STDERR( "[POINT CLASSIFICATION] : ---------------------------------" )
			return point_classification::outside;
		}
	}
	
	template <typename CoordinatesType>
	bool
	equals(const CoordinatesType& c1, const CoordinatesType& c2)
	{
		return (c1[0] == c2[0]) && (c1[1] == c2[1]) && (c1[2] == c2[2]);
	}
} //geometry

} //hmira

#endif // __POINT_IN_POLYHEDRON_HPP__