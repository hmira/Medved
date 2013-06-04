#ifndef __POINT_IN_POLYHEDRON_HPP__
#define __POINT_IN_POLYHEDRON_HPP__

#include <iostream>
#include <unordered_map>
#include <tuple>

#include <hmira/range/faces.hpp>
#include <hmira/geometry/line_on_plane_projection.hpp>

namespace hmira
{

namespace geometry
{
	template <typename TPoint, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	int point_in_polyhedron(
		TMesh& m,
		TPoint tested_point,
		TPoint tested_direction)// = TMesh_Traits::Point(1, 1, 1))
	{
// 		typedef std::pair<float, TMesh_Traits::face_descriptor> DistFacePair;
		std::unordered_map<float, int> distances;

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
// 					OpenMesh::Vec3f proj1(7,7,7), proj2(7,7,7);
					
// 					auto reflected = hmira::geometry::line_on_plane_projection(
// 						a,
// 						b,
// 						c,
// 						TMesh_Traits::get_face_normal(m, face),
// 						intersect,
// 						tested_direction,
// 						proj1,
// 						proj2);
					
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
					
// 					std::cout << "--" << dot_prod << " ref: " << reflected << " projected: " << proj1 << " : " << proj2 << std::endl;
					
					distances[dist] |= entry_flag;
					
// 					auto iter = distances.find(dist);
// 					std::cout << "flag: " << entry_flag << (distances[dist] | entry_flag) << std::endl;
// 					if (iter == distances.end())
// 					{
// 						distances[dist] = entry_flag;
// 					}
// 					else
// 					{
// 						iter->second |= entry_flag;
// 					}
				}
			});
		for (auto a:distances)
		{
			std::cout << a.second << std::endl;
		}
		
		auto boundary = false;
		int number_of_intersects = 0;
		
		for (auto a:distances)
		{
			if ( a.first <= 0.0001f &&  a.first >= -0.0001f )
				boundary = true;

			if ( a.second == 1 || a.second == 2 )
				number_of_intersects++;
		}
		
		if ( number_of_intersects & 1 )
		{
			std::cout << "inside" << std::endl;
			return 1; //inside
		}
		else
		{
			std::cout << "outside" << std::endl;
			return 0; //outside
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