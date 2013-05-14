#include "traits.h" //ak to nie je uplne hore, robi to nesmierne peklo :-)

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <algorithms/marching_cubes.hpp>

#include <meshes/OpenMeshX.hpp>

#include <grids/ScalarGridT.hh>
#include <grids/ScalarGrid_traits.h>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char **argv)
{	
	po::options_description desc("Allowed parameters");
	desc.add_options()
	("help,h","produce help message")
	("rasterize,r", po::value<std::string>()->default_value("full"), "type of rasterization [full|faces|edges]")
	("input-header-file,t", po::value<std::string>(), "input grid header file")
	("input-dump-file,i", po::value<std::string>(), "input grid dump file")
	("output-file,o", po::value<std::string>()->default_value("output.obj"), "output mesh .obj file")
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
	
	std::string input_header_filename;
	std::string input_dump_filename;
	
	int x,y,z;
	float x_size, y_size, z_size;

	if (vm.count("help"))
	{
		std::cerr << desc << std::endl;
		return 0;
	}
	
	if (vm.count("input-header-file"))
	{
		input_header_filename = vm["input-header-file"].as<std::string>();
	}
	else
	{
		std::cerr << "provide and input file\n" << desc << std::endl;
		return -1;
	}
	
	if (vm.count("input-dump-file"))
	{
		input_dump_filename = vm["input-dump-file"].as<std::string>();
	}
	else
	{
		std::cerr << "provide and input file\n" << desc << std::endl;
		return -1;
	}
	
	auto rasterization = vm["rasterize"].as<std::string>();
	if (rasterization != "full" && rasterization != "faces" && rasterization != "edges")
	{
		std::cerr << "unknown type of rasterization: \"" << rasterization << "\"\n" << desc << std::endl;
		return -1;
	}
	
	auto output_filename = vm["output-file"].as<std::string>();

	x = vm["x-resolution"].as<int>();
	y = vm["y-resolution"].as<int>();
	z = vm["z-resolution"].as<int>();
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

	OpenMeshExtended output_mesh;
	auto mc = MarchingCubes<IsoEx::ScalarGridT<float>, OpenMeshExtended, ScalarGrid_traits<float, IsoEx::ScalarGridT>>(sg, output_mesh);
	mc.process();
	
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