#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "marching_cubes_table.h"

template <
	typename TGrid,
	typename TMesh,
	typename TGrid_Traits,// = grid_traits<TGrid>, 
	typename TMesh_Traits = mesh_traits<TMesh>
	>
class MarchingCubes
{

public:
	typedef TGrid						Grid;

	typedef typename TGrid_Traits::Point_descriptor	Point_descriptor;
	typedef typename TGrid_Traits::Cube_descriptor		Cube_descriptor;
	typedef typename TGrid_Traits::Cube_iterator		Cube_iterator;

	typedef typename TMesh_Traits::Point			Mesh_Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;

	const Grid&	grid_;
	TMesh&            mesh_;

	std::vector<std::tuple<Point_descriptor, Point_descriptor, Vertex_descriptor>> vertices_;

	MarchingCubes( const Grid& _grid, TMesh& _mesh)
		: grid_( _grid ),
		mesh_( _mesh )
	{
		for ( auto cube : grid_ )
			process_cube( cube );
	}

	int marching_cubes(Grid& g, TMesh& m)
	{
		for ( auto cube : grid_ )
			process_cube( cube );
	}

	Vertex_descriptor
	add_vertex( Point_descriptor _p0, Point_descriptor _p1)
	{
		const Mesh_Point&  p0( TGrid_Traits::get_coords(grid_, _p0 ));
		const Mesh_Point&  p1( TGrid_Traits::get_coords(grid_, _p1 ));

		Vertex_descriptor vh;

		Mesh_Point p = Mesh_Point((p0 + ( p1-p0 )));

		if (_p1 < _p0) std::swap(_p0, _p1);

		for (auto itr : vertices_)
			if ( 	
				(std::get<0>(itr) == _p0) &&
				(std::get<1>(itr) == _p1) &&
				std::get<2>(itr).is_valid()
			)
				return std::get<2>(itr);

		vh = TMesh_Traits::create_vertex(p, mesh_);
		if (p0 < p1)
			vertices_.push_back( std::make_tuple(_p0, _p1, vh) );
		else
			vertices_.push_back( std::make_tuple(_p1, _p0, vh) );

		return vh;
	}

	void
	process_cube( Cube_descriptor _cidx )
	{

		Point_descriptor	corner[8];
		Vertex_descriptor	samples[12];
		unsigned int		cubetype( 0 );
		unsigned int		i;


		// get point indices of corner vertices
		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, _cidx, i );

		// determine cube type
		for ( i=0; i<8; ++i )
			if ( !TGrid_Traits::is_inside( grid_, corner[i] ))
				cubetype |= ( 1<<i );

		// trivial reject ?
		if ( cubetype == 0 || cubetype == 255 )
			return;

		// compute samples on cube's edges
		if ( edgeTable[cubetype]&1 )    samples[0]  = add_vertex( corner[0], corner[1] );
		if ( edgeTable[cubetype]&2 )    samples[1]  = add_vertex( corner[1], corner[2] );
		if ( edgeTable[cubetype]&4 )    samples[2]  = add_vertex( corner[3], corner[2] );
		if ( edgeTable[cubetype]&8 )    samples[3]  = add_vertex( corner[0], corner[3] );
		if ( edgeTable[cubetype]&16 )   samples[4]  = add_vertex( corner[4], corner[5] );
		if ( edgeTable[cubetype]&32 )   samples[5]  = add_vertex( corner[5], corner[6] );
		if ( edgeTable[cubetype]&64 )   samples[6]  = add_vertex( corner[7], corner[6] );
		if ( edgeTable[cubetype]&128 )  samples[7]  = add_vertex( corner[4], corner[7] );
		if ( edgeTable[cubetype]&256 )  samples[8]  = add_vertex( corner[0], corner[4] );
		if ( edgeTable[cubetype]&512 )  samples[9]  = add_vertex( corner[1], corner[5] );
		if ( edgeTable[cubetype]&1024 ) samples[10] = add_vertex( corner[2], corner[6] );
		if ( edgeTable[cubetype]&2048 ) samples[11] = add_vertex( corner[3], corner[7] );


		// connect samples by triangles
		for ( i=0; triTable[cubetype][0][i] != -1; i+=3 )
		TMesh_Traits::create_face(
			samples[triTable[cubetype][0][i  ]],
			samples[triTable[cubetype][0][i+1]],
			samples[triTable[cubetype][0][i+2]],
			mesh_ );
	}
};


