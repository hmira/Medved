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
	("help","produce help message")
	("input-file,i", po::value<std::string>(), "input file")
	("output-file,o", po::value<std::string>(), "output file")
	("x-resolution,x", po::value<int>(), "x resolution");
	
	po::positional_options_description p;
	p.add("input-file,i", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
		options(desc).positional(p).run(), vm);
	po::notify(vm);
	
	std::string input_filename;
	std::string output_filename;
	
	std::cout<<"ty si taký kokot"<<std::endl;
	
	if (vm.count("input-file"))
	{
		input_filename = vm["input-file"].as<std::string>();
	}
	else
	{
		return 0;
	}
	
	if (vm.count("output-file"))
	{
		output_filename = vm["output-file"].as<std::string>();
	}
	else
	{
		return 0;
	}
	
	if (vm.count("x-resolution"))
	{
		std::cout << "x:" << vm["x-resolution"].as<int>() << std::endl;
	}
	else
	{
		std::cout << "dopiči" << std::endl;
	}
	
	OpenMeshExtended input_mesh, output_mesh;
	if (!OpenMesh::IO::read_mesh(input_mesh, input_filename))
	{
		std::cerr << "error reading file:" << input_filename << std::endl;
		return 1;
	}


	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
OpenMesh::VectorT<float, 3>( 3, 0, 0 ),
OpenMesh::VectorT<float, 3>( 0, 3, 0 ),
OpenMesh::VectorT<float, 3>( 0, 0, 3 ),
30,
30,
30);
	
	auto vx = Voxelize<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, input_mesh);	
 	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, output_mesh);
// 	fill_holes<OpenMeshExtended, advanced_mesh_traits<OpenMeshExtended>>(mesh2);



	
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