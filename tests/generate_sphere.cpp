#include <boost/mpl/if.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/int.hpp>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include <iostream>
#include <deque>
#include <hmira/meshes/winged_edge.hpp>

#include <hmira/algorithms/compute_components.hpp>
#include <hmira/algorithms/flip_normals.hpp>
#include <hmira/algorithms/triangulate.hpp>

#include <hmira/algorithms/marching_cubes.hpp>
#include <hmira/algorithms/voxelize.hpp>
#include <hmira/algorithms/fill_holes.hpp>

#include <hmira/meshes/OpenMeshX.hpp>
#include "traits.h"

#include <hmira/grids/ScalarGridT.hh>
#include <hmira/grids/ScalarGrid_traits.h>
#include <hmira/grids/ImplicitSphere.hh>

#include <math.h>
#include <hmira/geometry/ray_face_intersection.hpp>

template <typename T>
constexpr auto has_x_method(T& t) -> decltype(t.add_vertex(), OpenMeshExtended::VertexHandle())
{
	return true;
}

typedef winged_edge_mesh<triangleMesh> my_mesh;
typedef winged_edge_mesh_traits<triangleMesh> my_mesh_traits;

int main(int argc, char **argv)
{
	po::options_description desc("Allowed parameters");
	desc.add_options()
	("help,h","produce help message")
	("output-raw-file,o", po::value<std::string>()->default_value("output.dump"), "output grid raw file")
	("output-header-file,t", po::value<std::string>()->default_value("output.hdr"), "output grid header file")
	("x-resolution,x", po::value<int>()->default_value(30), "x resolution")
	("y-resolution,y", po::value<int>()->default_value(30), "y resolution")
	("z-resolution,z", po::value<int>()->default_value(30), "z resolution")
	("x-size,X", po::value<float>()->default_value(3.f), "X size")
	("y-size,Y", po::value<float>()->default_value(3.f), "Y size")
	("z-size,Z", po::value<float>()->default_value(3.f), "Z size");
	
	po::positional_options_description p;
	p.add("input-file,i", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
		options(desc).positional(p).run(), vm);
	po::notify(vm);
	
	std::string input_filename;
	
	int x,y,z;
	float x_size, y_size, z_size;

	if (vm.count("help"))
	{
		std::cerr << desc << std::endl;
		return 0;
	}
	
	auto output_raw_filename = vm["output-raw-file"].as<std::string>();
	auto output_header_filename = vm["output-header-file"].as<std::string>();

	x = vm["x-resolution"].as<int>();
	y = vm["y-resolution"].as<int>();
	z = vm["z-resolution"].as<int>();
	x_size = vm["x-size"].as<float>();
	y_size = vm["y-size"].as<float>();
	z_size = vm["z-size"].as<float>();
	
		
	OpenMeshExtended input_mesh, output_mesh;
	if (!OpenMesh::IO::read_mesh(input_mesh, input_filename))
	{
		std::cerr << "error reading file:" << input_filename << std::endl;
		return 1;
	}
	
	std::cerr << 	"[GRID] : Buiding empty ScalarGridT<float>\nresolution:\n" <<
			"x: " << x << " cubes\n" <<
			"y: " << y << " cubes\n" <<
			"z: " << z << " cubes" << std::endl;
	std::cerr << "Grid size:\n" <<
			"x: " << x_size << " cubes\n" <<
			"y: " << y_size << " cubes\n" <<
			"z: " << z_size << " cubes" << std::endl;

	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
		OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
		OpenMesh::VectorT<float, 3>( x_size, 0, 0 ),
		OpenMesh::VectorT<float, 3>( 0, y_size, 0 ),
		OpenMesh::VectorT<float, 3>( 0, 0, z_size ),
		x,
		y,
		z);
	
	
	auto origin = OpenMesh::Vec3f(0, 0.3, 0.3);
	auto direction = OpenMesh::Vec3f(1,0,0);


	
	IsoEx::ImplicitSphere ims = IsoEx::ImplicitSphere(OpenMesh::Vec3f(1.5f,1.5f,1.5f), 1.5f);
	sg.sample(ims);

	std::cerr << "[GRID] : write to header file \"" << output_header_filename << "\"" << std::endl;
	ScalarGrid_traits<float, IsoEx::ScalarGridT>::write_header(output_header_filename, sg);
	std::cerr << "[GRID] : write to raw file \"" << output_raw_filename << "\"" << std::endl;
	ScalarGrid_traits<float, IsoEx::ScalarGridT>::write_dump(output_raw_filename, sg);
}