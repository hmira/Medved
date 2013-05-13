#include <boost/mpl/if.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/int.hpp>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include <iostream>
#include <deque>
#include <string>
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

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv)
{	
	po::options_description desc("Allowed parameters");
	desc.add_options()
	("help,h","produce help message")
	("fill-holes,f","fill holes on the resulting mesh")
	("input-file,i", po::value<std::string>(), "input file")
	("output-file,o", po::value<std::string>()->default_value("output.obj"), "output file")
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
	
	if (vm.count("input-file"))
	{
		input_filename = vm["input-file"].as<std::string>();
	}
	else
	{
		std::cerr << "provide and input file\n" << desc << std::endl;
		return -1;
	}
	
	auto fill = (vm.count("fill-holes"));
	
	auto output_filename = vm["output-file"].as<std::string>();

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
	
	std::cerr << 	"[GRID]Buiding empty ScalarGridT<float>\nresolution:\n" <<
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
	
	auto vx = Voxelize<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, input_mesh);	
 	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, output_mesh);

	if (fill)
	{
		fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(output_mesh);
	}


	
	if (!OpenMesh::IO::write_mesh(output_mesh, output_filename)) 
	{
		std::cerr << "write error\n";
		exit(1);
	}
	else
	{
		std::cerr << "object: " << output_filename << " written" <<std::endl;
	}


	return 0;
}