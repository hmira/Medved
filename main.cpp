/*
 * main.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: hmirap
 */

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>


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


	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
OpenMesh::VectorT<float, 3>( 1, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 1, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 1 ),
50,
50,
50 
);


	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(0.5f,0.5f,0.5f), 0.5f);

	sg.sample(ims);

	MarchingCubes<OpenMeshExtended> mc = MarchingCubes<OpenMeshExtended>(sg, mesh);
	

	std::cout << "norma: c++0x" << std::endl;

	std::cout << sg.n_cubes() << " aa " << sg.n_points() << std::endl;



	auto my_pair_ww = OpenMeshXTraits::get_all_vertices(mesh);
	auto my_pair_ff = OpenMeshXTraits::get_all_faces(mesh);
	int pocitam = 0;
	int pocitam_faces = 0;
	for (auto i = my_pair_ww.first; i != my_pair_ww.second; ++i) {
//		std::cout << "aleluja" << ", ";
		pocitam++;
	}
	for (auto i = my_pair_ff.first; i != my_pair_ff.second; ++i) {
		pocitam_faces++;
	}
	std::cout << std::endl;
	std::cout << "pocitam vertices: " << pocitam << std::endl;
	std::cout << "pocitam faces: " << pocitam_faces << std::endl;

	if (!OpenMesh::IO::write_mesh(mesh, "my_mesh.obj")) 
	{
		std::cerr << "write error\n";
		exit(1);
	}

	//M4D::Imaging::AImage::Ptr image = M4D::Imaging::ImageFactory::LoadDumpedImage( path );

	return 0;
}
