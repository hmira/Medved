#include "traits.h" //ak to nie je uplne hore, robi to nesmierne peklo :-)

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <hmira/algorithms/marching_cubes.hpp>
#include <hmira/algorithms/fill_holes.hpp>

#include <hmira/meshes/OpenMeshX.hpp>

#include <hmira/grids/ScalarGridT.hh>
#include <hmira/grids/ScalarGrid_traits.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv)
{	
	po::options_description desc("Allowed parameters");
	desc.add_options()
	("help,h","produce help message")
	("input-header-file,t", po::value<std::string>(), "input grid header file")
	("input-dump-file,i", po::value<std::string>(), "input grid dump file")
	("output-file,o", po::value<std::string>()->default_value("output.obj"), "output mesh .obj file")
	("x-size,X", po::value<float>()->default_value(3.f), "X size")
	("y-size,Y", po::value<float>()->default_value(3.f), "Y size")
	("z-size,Z", po::value<float>()->default_value(3.f), "Z size");
	
	po::positional_options_description p;
	p.add("input-file,i", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
		options(desc).positional(p).run(), vm);
	po::notify(vm);
	
	std::string input_header_filename;
	std::string input_dump_filename;
	
	int x,y,z;
	float x_size, y_size, z_size;

	if (vm.count("help"))
	{
		std::cerr << desc << std::endl;
		return 0;
	}
	
	if ((vm.count("input-header-file")) && (vm.count("input-dump-file")))
	{
		input_header_filename = vm["input-header-file"].as<std::string>();
		input_dump_filename = vm["input-dump-file"].as<std::string>();
	}
	else
	{
		std::cerr << "provide and input file\n" << desc << std::endl;
		return -1;
	}
		
	auto output_filename = vm["output-file"].as<std::string>();

	x_size = vm["x-size"].as<float>();
	y_size = vm["y-size"].as<float>();
	z_size = vm["z-size"].as<float>();
	
	ScalarGrid_traits<float, IsoEx::ScalarGridT>::read_header(input_header_filename, x, y, z);
	
	IsoEx::ScalarGridT<float> sg = IsoEx::ScalarGridT<float>(
	OpenMesh::VectorT<float, 3>( 0, 0, 0 ),
	OpenMesh::VectorT<float, 3>( x_size, 0, 0 ),
	OpenMesh::VectorT<float, 3>( 0, y_size, 0 ),
	OpenMesh::VectorT<float, 3>( 0, 0, z_size ),
	x + 1,
	y + 1,
	z + 1);
	
	ScalarGrid_traits<float, IsoEx::ScalarGridT>::read_dump(input_dump_filename, sg);
	std::cerr << "[GRID] : file " << input_dump_filename << " loaded." << std::endl;

	OpenMeshExtended output_mesh;
	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, output_mesh);
	mc.parallel_process();
	
	if (!OpenMesh::IO::write_mesh(output_mesh, output_filename)) 
	{
		std::cerr << "write error\n";
		exit(1);
	}
	else
	{
		std::cerr << "[MESH] : object: " << output_filename << " written" <<std::endl;
	}


}