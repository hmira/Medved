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
#include "algorithms/fill_holes.h"

#include "meshes/OpenMeshX.h"
#include "traits.h"

#include <grids/ScalarGridT.hh>
#include <grids/ScalarGrid_traits.h>
#include <grids/ImplicitSphere.hh>

typedef winged_edge_mesh<triangleMesh> my_mesh;
typedef winged_edge_mesh_traits<triangleMesh> my_mesh_traits;

enum {A,B,C,D,N};

int main(int argc, char **argv)
{
	OpenMeshExtended mesh;
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
 // mesh.add_face(face_vhandles);

	std::cout << vhandle[0] << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.halfedge_handle(vhandle[0])) << std::endl;
	std::cout << mesh.to_vertex_handle(mesh.opposite_halfedge_handle(mesh.halfedge_handle(vhandle[0]))) << std::endl;

//	advanced_mesh_traits<OpenMeshExtended>::truncate(mesh, vhandle[6], 0.2f);
//	advanced_mesh_traits<OpenMeshExtended>::truncate(mesh, (OpenMeshExtended::VertexHandle)8);
*/


	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
OpenMesh::VectorT<float, 3>( 1, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 1, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 1 ),
25,
25,
25 
);


	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(0.5f,0.5f,0.5f), 0.6f);

	sg.sample(ims);

	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, mesh);
	mesh.update_face_normals();	

	fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	triangulate<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh);
	
	std::cout << "norma: c++0x" << std::endl;
	
	if (!OpenMesh::IO::write_mesh(mesh, "my_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}

	//M4D::Imaging::AImage::Ptr image = M4D::Imaging::ImageFactory::LoadDumpedImage( path );

	return 0;
}
