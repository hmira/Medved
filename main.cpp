/*
 * main.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: hmirap
 */

#include <boost/mpl/if.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/int.hpp>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
//#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

#include <iostream>
#include <deque>
#include "meshes/winged_edge.hpp"

#include "algorithms/compute_components.hpp"
#include "algorithms/flip_normals.hpp"
#include "algorithms/triangulate.hpp"

#include "algorithms/marching_cubes.hpp"
#include "algorithms/voxelize.hpp"
#include "algorithms/fill_holes.hpp"

#include "meshes/OpenMeshX.hpp"
#include "traits.h"

#include <grids/ScalarGridT.hh>
#include <grids/ScalarGrid_traits.h>
#include <grids/ImplicitSphere.hh>

#include <math.h>
#include <geometry/ray_face_intersection.hpp>

template <typename T>
constexpr auto has_x_method(T& t) -> decltype(t.add_vertex(), OpenMeshExtended::VertexHandle())
{
	return true;
}

typedef winged_edge_mesh<triangleMesh> my_mesh;
typedef winged_edge_mesh_traits<triangleMesh> my_mesh_traits;

int main(int argc, char **argv)
{
	auto origin = OpenMesh::Vec3f(0, 0.3, 0.3);
	auto direction = OpenMesh::Vec3f(1,0,0);

	auto v1 = OpenMesh::Vec3f(1,0,0);
	auto v2 = OpenMesh::Vec3f(1,2,0);
	auto v3 = OpenMesh::Vec3f(1,0,1);	
	
	OpenMesh::Vec3f intersection(0,0,0);
	float distance = 0.f;
	
	auto resultt = line_face_intersection(
		origin,
		direction,
		v1,
		v2,
		v3,
		intersection,
		distance
	);
	/*
	std::cout << resultt << std::endl;
	std::cout << intersection << std::endl;
	return 0;*/
	
	OpenMeshExtended mesh, mesh2;

	std::cout << std::boolalpha << std::is_member_function_pointer<decltype(&OpenMeshExtended::add_vertex)>::value << std::endl;

	auto point1 = OpenMeshExtended::Point(-1,-1,1);
	auto point2 = OpenMeshExtended::Point(-2,1,0.2);
	auto res = cross(point1, point2);
	std::cout << res << std::endl;

	OpenMeshExtended::VertexHandle vhandle[4];
	vhandle[0] = mesh.add_vertex(OpenMeshExtended::Point(-1.1, -1.1,  -1.1));
	vhandle[1] = mesh.add_vertex(OpenMeshExtended::Point( 1.1, -1.1,  -1.1));
	vhandle[2] = mesh.add_vertex(OpenMeshExtended::Point( 0.1, -0.5, 1.1));
	vhandle[3] = mesh.add_vertex(OpenMeshExtended::Point( 0.1,  1.1,   0.1));
	auto f0 = mesh.add_face(vhandle[3], vhandle[1], vhandle[0]);
	auto f1 = mesh.add_face(vhandle[3], vhandle[2], vhandle[1]);
	auto f2 = mesh.add_face(vhandle[3], vhandle[0], vhandle[2]);
	auto f3 = mesh.add_face(vhandle[0], vhandle[1], vhandle[2]);
	
	std::cout << "normal: " << mesh_traits<OpenMeshExtended>::get_face_normal(mesh, f0) << std::endl;
	std::cout << "normal: " << mesh_traits<OpenMeshExtended>::get_face_normal(mesh, f1) << std::endl;
	std::cout << "normal: " << mesh_traits<OpenMeshExtended>::get_face_normal(mesh, f2) << std::endl;
	std::cout << "normal: " << mesh_traits<OpenMeshExtended>::get_face_normal(mesh, f3) << std::endl;
	
//   OpenMeshExtended::VertexHandle vhandle[8];
//   vhandle[0] = mesh.add_vertex(OpenMeshExtended::Point(-1, -1,  1));
//   vhandle[1] = mesh.add_vertex(OpenMeshExtended::Point( 1, -1,  1));
//   vhandle[2] = mesh.add_vertex(OpenMeshExtended::Point( 1,  1,  1));
//   vhandle[3] = mesh.add_vertex(OpenMeshExtended::Point(-1,  1,  1));
//   vhandle[4] = mesh.add_vertex(OpenMeshExtended::Point(-1, -1, -1));
//   vhandle[5] = mesh.add_vertex(OpenMeshExtended::Point( 1, -1, -1));
//   vhandle[6] = mesh.add_vertex(OpenMeshExtended::Point( 1,  1, -1));
//   vhandle[7] = mesh.add_vertex(OpenMeshExtended::Point(-1,  1, -1));
//   // generate (quadrilateral) faces
//   std::vector<OpenMeshExtended::VertexHandle>  face_vhandles;
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[0]);
//   face_vhandles.push_back(vhandle[1]);
//   face_vhandles.push_back(vhandle[2]);
//   face_vhandles.push_back(vhandle[3]);
//   mesh.add_face(face_vhandles);
//  
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[7]);
//   face_vhandles.push_back(vhandle[6]);
//   face_vhandles.push_back(vhandle[5]);
//   face_vhandles.push_back(vhandle[4]);
//   mesh.add_face(face_vhandles);
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[1]);
//   face_vhandles.push_back(vhandle[0]);
//   face_vhandles.push_back(vhandle[4]);
//   face_vhandles.push_back(vhandle[5]);
//   mesh.add_face(face_vhandles);
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[2]);
//   face_vhandles.push_back(vhandle[1]);
//   face_vhandles.push_back(vhandle[5]);
//   face_vhandles.push_back(vhandle[6]);
//   mesh.add_face(face_vhandles);
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[3]);
//   face_vhandles.push_back(vhandle[2]);
//   face_vhandles.push_back(vhandle[6]);
//   face_vhandles.push_back(vhandle[7]);
//   mesh.add_face(face_vhandles);
//   face_vhandles.clear();
//   face_vhandles.push_back(vhandle[0]);
//   face_vhandles.push_back(vhandle[3]);
//   face_vhandles.push_back(vhandle[7]);
//   face_vhandles.push_back(vhandle[4]);
	
/*
  OpenMeshExtended::VertexHandle vhandle[8];
  vhandle[0] = mesh.add_vertex(OpenMeshExtended::Point(-1, -1,  1));
  vhandle[1] = mesh.add_vertex(OpenMeshExtended::Point( 1, -1,  1));
  vhandle[2] = mesh.add_vertex(OpenMeshExtended::Point( 1,  1,  1));
  vhandle[3] = mesh.add_vertex(OpenMeshExtended::Point(-1,  1,  1));
  vhandle[4] = mesh.add_vertex(OpenMeshExtended::Point(-1, -1, -1));
  vhandle[5] = mesh.add_vertex(OpenMeshExtended::Point( 1, -1, -1));
  vhandle[6] = mesh.add_vertex(OpenMeshExtended::Point( 1,  1, -1));
  vhandle[7] = mesh.add_vertex(OpenMeshExtended::Point(-1,  1, -1));
  // generate (quadrilateral) faces
  std::vector<OpenMeshExtended::VertexHandle>  face_vhandles;
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[0]);
  face_vhandles.push_back(vhandle[1]);
  face_vhandles.push_back(vhandle[2]);
  face_vhandles.push_back(vhandle[3]);
  mesh.add_face(face_vhandles);
 
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[7]);
  face_vhandles.push_back(vhandle[6]);
  face_vhandles.push_back(vhandle[5]);
  face_vhandles.push_back(vhandle[4]);
  mesh.add_face(face_vhandles);
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[1]);
  face_vhandles.push_back(vhandle[0]);
  face_vhandles.push_back(vhandle[4]);
  face_vhandles.push_back(vhandle[5]);
  mesh.add_face(face_vhandles);
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[2]);
  face_vhandles.push_back(vhandle[1]);
  face_vhandles.push_back(vhandle[5]);
  face_vhandles.push_back(vhandle[6]);
  mesh.add_face(face_vhandles);
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[3]);
  face_vhandles.push_back(vhandle[2]);
  face_vhandles.push_back(vhandle[6]);
  face_vhandles.push_back(vhandle[7]);
  mesh.add_face(face_vhandles);
  face_vhandles.clear();
  face_vhandles.push_back(vhandle[0]);
  face_vhandles.push_back(vhandle[3]);
  face_vhandles.push_back(vhandle[7]);
  face_vhandles.push_back(vhandle[4]);
>>>>>>> 30d3b21ed3781ae60c4b815f9df849fd6318078d
 // mesh.add_face(face_vhandles);
/*
	std::cout << vhandle[0] << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.halfedge_handle(vhandle[0])) << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.opposite_halfedge_handle(mesh.halfedge_handle(vhandle[0]))) << std::endl;

	auto scale = advanced_mesh_traits<OpenMeshExtended>::bevel(mesh, (advanced_mesh_traits<OpenMeshExtended>::edge_descriptor)0, 0.5f);
//	advanced_mesh_traits<OpenMeshExtended>::truncate(mesh, vhandle[5], 0.5f);

//	scale.rescale(0.8f);
	scale.rescale_by_unit(1.0f);;
	scale.rescale_by_unit(-1.0f);;
	scale.rescale_by_unit(-0.3f);;

//	scale.rescale_by_unit(0.5f);
//	scale.rescale_by_unit(0.5f);
	
*/
	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
//OpenMesh::VectorT<float, 3>( -2.2, -2.2, -2.2 ),
OpenMesh::VectorT<float, 3>( -1.5, -1.5, -1.5 ),
OpenMesh::VectorT<float, 3>( 3, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 3, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 3 ),
30,
30,
30);
// */


	
//	auto vx = Voxelize<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, mesh);


	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(0.f,0.f,0.f), 1.4f);

	sg.sample(ims);
	
	auto firstpair = ScalarGrid_traits<float, IsoEx::ScalarGridT>::get_bounds(sg, 2, 0);
	auto secondpair = ScalarGrid_traits<float, IsoEx::ScalarGridT>::get_bounds(sg, 2, 1);

	std::cout << "1pair: " << firstpair.first << " & " << firstpair.second << std::endl;
	std::cout << "2pair: " << secondpair.first << " & " << secondpair.second  << std::endl;
	
	
	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, mesh2);
//	mesh.update_face_normals();	

	
/*	fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);

	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, mesh);
	mesh.update_face_normals();	

	fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	triangulate<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	
	advanced_mesh_traits<OpenMeshExtended>::flip_face_normal(
		mesh,
		(advanced_mesh_traits<OpenMeshExtended>::face_descriptor)0
								);
*/

// 	auto all_faces = mesh_traits<OpenMeshExtended>::get_all_faces(mesh2);
// 	for (auto f = all_faces.first; f != all_faces.second; ++f)
// 	{
// 		std::cout << "OO" << std::endl;
// 
// 		auto my_face = *f;
// 		auto surr_v = mesh_traits<OpenMeshExtended>::get_surrounding_vertices(mesh, my_face);
// 
// 		for (auto it_sv = surr_v.first; it_sv != surr_v.second; ++it_sv)
// 		{
// 		}
// 	}

	
	std::cout << "norma: c++0x" << std::endl;
	
	if (!OpenMesh::IO::write_mesh(mesh, "my_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}
	
	if (!OpenMesh::IO::write_mesh(mesh2, "snd_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}

	//M4D::Imaging::AImage::Ptr image = M4D::Imaging::ImageFactory::LoadDumpedImage( path );

	return 0;
}
