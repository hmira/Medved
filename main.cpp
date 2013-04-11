/*
 * main.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: hmirap
 */

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>

//#include "MedV4D/Imaging/ImageTools.h"


typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

#include <iostream>
#include <deque>
#include <boost/graph/adjacency_matrix.hpp>
#include "meshes/winged_edge.h"

#include "algorithms/compute_components.h"
#include "algorithms/flip_normals.h"
#include "algorithms/triangulate.h"

#include "algorithms/marching_cubes.h"

#include "meshes/OpenMeshX.h"
#include "traits.h"

#include <grids/ScalarGridT.hh>
#include <grids/ImplicitSphere.hh>

typedef winged_edge_mesh<triangleMesh> my_mesh;
typedef winged_edge_mesh_traits<triangleMesh> my_mesh_traits;

enum {A,B,C,D,N};

int main(int argc, char **argv)
{
	OpenMeshExtended mesh;

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
  mesh.add_face(face_vhandles);

	std::cout << vhandle[0] << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.halfedge_handle(vhandle[0])) << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.opposite_halfedge_handle(mesh.halfedge_handle(vhandle[0]))) << std::endl;

	advanced_mesh_traits<OpenMeshExtended>::truncate(mesh, vhandle[0]);

/*
	auto e_itr   = mesh.edges_begin();
	auto e_end   = mesh.edges_end();
	for ( ; e_itr != e_end; ++e_itr)
		advanced_mesh_traits<OpenMeshExtended>::split_edge(mesh, e_itr.handle());


	e_itr   = mesh.edges_begin();
	e_end   = mesh.edges_end();
	for ( ; e_itr != e_end; ++e_itr)
		advanced_mesh_traits<OpenMeshExtended>::split_edge(mesh, e_itr.handle());


	for ( ; e_itr != e_end; ++e_itr)
		std::cout << mesh.is_boundary(e_itr.handle());
*/
	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
OpenMesh::VectorT<float, 3>( 1, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 1, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 1 ),
5,
5,
5 
);


	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(0.5f,0.5f,0.5f), 0.6f);

	sg.sample(ims);

//	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended>(sg, mesh);
//	mesh.update_face_normals();	

	std::cout << "norma: c++0x" << std::endl;

	std::cout << sg.n_cubes() << " aa " << sg.n_points() << std::endl;



	auto my_pair_ww = OpenMeshXTraits::get_all_vertices(mesh);
	auto my_pair_ff = OpenMeshXTraits::get_all_faces(mesh);
	int pocitam = 0;
	int pocitam_faces = 0;
	for (auto i = my_pair_ww.first; i != my_pair_ww.second; ++i) {
		pocitam++;
	}
	for (auto i = my_pair_ff.first; i != my_pair_ff.second; ++i) {
		mesh.set_color(i.handle(), OpenMeshExtended::Color(1,0,0));
		pocitam_faces++;
		int fv_count = 0;
		auto fv_pair = OpenMeshXTraits::get_surrounding_vertices
			(mesh, 
			i.handle());
		for (auto j = fv_pair.first; j!=fv_pair.second; ++j)
			fv_count++;

		std::cout << "vrcholov: " << fv_count << std::endl;
	}

	std::cout << std::endl;
	std::cout << "pocitam vertices: " << pocitam << std::endl;
	std::cout << "pocitam faces: " << pocitam_faces << std::endl;

	OpenMesh::Subdivider::Uniform::CatmullClarkT<OpenMeshExtended> catmull;
	

 	catmull.attach(mesh);
//	catmull( 2 );
	catmull.detach();


//	triangulate<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);

	//advanced_mesh_traits<OpenMeshExtended>::flip_face_normal_t<int>();

	//mesh.update_face_normals();
	//flip_normals<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);	

	/*auto ff_face_pair = mesh_traits<OpenMeshExtended>::get_all_faces(mesh);
	for (auto i = ff_face_pair.first; i != ff_face_pair.second; ++i)
	{
		auto old_face = *i;
		OpenMeshExtended::Normal n = mesh.calc_face_normal(old_face);
		std::cout << "norm:" << n << std::endl;
		std::cout << "n_normal:" << mesh.normal(old_face) << std::endl;
	} */

	if (!OpenMesh::IO::write_mesh(mesh, "my_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}

	//M4D::Imaging::AImage::Ptr image = M4D::Imaging::ImageFactory::LoadDumpedImage( path );

	return 0;
}
