#include "traits.h" //ak to nie je uplne hore, robi to nesmierne peklo :-)

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <iostream>
#include <string>

#include <hmira/meshes/OpenMeshX.hpp>
#include <hmira/geometry/point_in_polyhedron.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/**
 * \brief the main method of Point in Polyhedron test program
 * 
 * The allowed parameters from command line are shown
 * after option --help
 * 
 * Polyhedron is represented by a standard polygonal mesh
 * and point is represented by a triplet of coordinates in 3D
 * 
 */
int main(int argc, char **argv)
{
	float x, y, z;
	std::string input_filename;
	po::options_description desc("Allowed parameters");
	desc.add_options()
	("help,h","produce help message")
	("input-file,i", po::value<std::string>(), "input mesh file")
	("x-coordinate,X", po::value<float>(&x)->default_value(0.f), "X coordinate")
	("y-coordinate,Y", po::value<float>(&y)->default_value(0.f), "Y coordinate")
	("z-coordinate,Z", po::value<float>(&z)->default_value(0.f), "Z coordinate");

	po::positional_options_description p;
	p.add("input-file,i", -1);

	try
	{
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);

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
	}
	catch(po::error& e)
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
		return 1;
	}
	
	auto point = OpenMeshExtended::Point(x,  y, z);
	OpenMeshExtended input_mesh;
	if (!OpenMesh::IO::read_mesh(input_mesh, input_filename))
	{
		std::cerr << "error reading file:" << input_filename << std::endl;
		return 1;
	}
	
	auto result = hmira::geometry::point_in_polyhedron(input_mesh, point);
	std::cout << "the point: " << point << " is classificated as: " << std::endl;

	if (result == hmira::geometry::point_classification::inside)
		std::cout << "INSIDE" << std::endl;
	else if (result == hmira::geometry::point_classification::outside)
		std::cout << "OUTSIDE" << std::endl;
	else if (result == hmira::geometry::point_classification::boundary)
		std::cout << "BOUNDARY";
}