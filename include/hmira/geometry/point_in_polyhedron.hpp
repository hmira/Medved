#ifndef __POINT_IN_POLYHEDRON_HPP__
#define __POINT_IN_POLYHEDRON_HPP__

#include <hmira/range/faces.hpp>

namespace hmira
{

namespace geometry
{
	template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	int point_in_polyhedron(TMesh& m)
	{
// 		range::for_each_face(m,
// 			[&]
// 			(TMesh_Traits::face_descriptor face)
// 			{
// 				auto fv_pair = TMesh_Traits::get_surrounding_vertices(m, face);
// 				TMesh_Traits::vertex_descriptor v[3];
// 				
// 				float dist = 0.f;
// 				OpenMesh::Vec3f intersect(0,0,0);
// 				
// 				auto fv = fv_pair.first;
// 				for(int i = 0; i < 3; i++ )
// 				{
// 					v[i] = *fv;
// 					++fv;
// 				}
// 				
// 				auto a = TMesh_Traits::get_coordinates(m, v[0]);
// 				auto b = TMesh_Traits::get_coordinates(m, v[1]);
// 				auto c = TMesh_Traits::get_coordinates(m, v[2]);
// 			
// 				auto result = line_face_intersection(
// 				OpenMesh::Vec3f(0, 0, -1),
// 				OpenMesh::Vec3f(0, 0, -1),
// 				a,
// 				b,
// 				c,
// 				intersect,
// 				dist);
// 				
// 				std::cout << "vec: " << a << std::endl;
// 				std::cout << "vec: " << b << std::endl;
// 				std::cout << "vec: " << c << std::endl;
// 				std::cout << "result: " << result << (intersect) << std::endl;
// 			});
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