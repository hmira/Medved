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
#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

#include <iostream>
#include <deque>
#include "meshes/winged_edge.hpp"

#include "algorithms/compute_components.hpp"
#include "algorithms/flip_normals.hpp"
#include "algorithms/triangulate.hpp"

#include "algorithms/marching_cubes.hpp"
#include "algorithms/fill_holes.hpp"

#include "meshes/OpenMeshX.hpp"
#include "traits.h"

#include <grids/ScalarGridT.hh>
#include <grids/ScalarGrid_traits.h>
#include <grids/ImplicitSphere.hh>

template <typename T>
constexpr auto has_x_method(T& t) -> decltype(t.add_vertex(), OpenMeshExtended::VertexHandle())
{
	return true;
}

typedef winged_edge_mesh<triangleMesh> my_mesh;
typedef winged_edge_mesh_traits<triangleMesh> my_mesh_traits;

int main(int argc, char **argv)
{

	OpenMeshExtended mesh;

	std::cout << std::boolalpha << std::is_member_function_pointer<decltype(&OpenMeshExtended::add_vertex)>::value << std::endl;


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
 // mesh.add_face(face_vhandles);

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

/*	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
OpenMesh::VectorT<float, 3>( 1, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 1, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 1 ),
10,
10,
10 
);


	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(0.5f,0.5f,0.5f), 0.5f);

	sg.sample(ims);

	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, mesh);
	mesh.update_face_normals();	
*/
	fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	triangulate<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	
	advanced_mesh_traits<OpenMeshExtended>::flip_face_normal(
		mesh,
		(advanced_mesh_traits<OpenMeshExtended>::face_descriptor)0
								);

/*
	auto all_faces = mesh_traits<OpenMeshExtended>::get_all_faces(mesh);
	for (auto f = all_faces.first; f != all_faces.second; ++f)
	{
		auto my_face = *f;
		auto surr_v = mesh_traits<OpenMeshExtended>::get_surrounding_vertices(mesh, my_face);

		for (auto it_sv = surr_v.first; it_sv != surr_v.second; ++it_sv)
		{
		}
	}
*/
	
	std::cout << "norma: c++0x" << std::endl;
	
	if (!OpenMesh::IO::write_mesh(mesh, "my_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}

	//M4D::Imaging::AImage::Ptr image = M4D::Imaging::ImageFactory::LoadDumpedImage( path );

	return 0;
}
