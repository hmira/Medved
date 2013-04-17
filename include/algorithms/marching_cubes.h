#include <grids/MCTablesIsoEx.hh>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>

template <typename TGrid, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
class MarchingCubes
{

public:
	typedef TGrid			Grid;
	typedef typename TGrid::PointIdx	PointIdx;
	typedef typename TGrid::CubeIdx		CubeIdx;
	typedef typename TGrid::CubeIterator	CubeIterator;

	typedef typename TMesh_Traits::Point			Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;

	const Grid&	grid_;
	TMesh&            mesh_;
	float            iso_;

	std::vector<std::tuple<PointIdx, PointIdx, Vertex_descriptor>> vertices_;

	MarchingCubes( const Grid& _grid, TMesh& _mesh, float _iso = 0.0f)
		: grid_( _grid ),
		mesh_( _mesh ),
		iso_( _iso )
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
	add_vertex( PointIdx _p0, PointIdx _p1)
	{
		const Point&  p0( grid_.point( _p0 ) );
		const Point&  p1( grid_.point( _p1 ) );

		Vertex_descriptor vh;		

		float s0 = grid_.scalar_distance( _p0 );
		float s1 = grid_.scalar_distance( _p1 );

		if ( fabs( s1-s0 ) > 0.0001f )
		{
			float t  = ( iso_-s0 ) / ( s1-s0 );
			Point p = Point((p0 + ( p1-p0 )*t));

			if (_p1 < _p0) std::swap(_p0, _p1);

			for (auto itr : vertices_)
				if ( 	
					((std::get<0>(itr) == _p0)) &&
					((std::get<1>(itr) == _p1)) &&
					std::get<2>(itr).is_valid()
				)
				{	
					return std::get<2>(itr);
				}


			//vh = mesh_.add_vertex( p );
			vh = TMesh_Traits::create_vertex(p, mesh_);
			if (p0 < p1)
				vertices_.push_back( std::make_tuple(_p0, _p1, vh) );
			else
				vertices_.push_back( std::make_tuple(_p1, _p0, vh) );

			return vh;
		}
		else
		{
			std::cout << s0 << ", " <<s1 <<"!!!!!!!!!!!!!!!!!" <<( s1 - s0 ) << std::endl;
			return vh;
		}
	}

	void
	process_cube( CubeIdx _cidx )
	{

	PointIdx		corner[8];
	Vertex_descriptor	samples[12];
	unsigned char		cubetype( 0 );
	unsigned int		i;


	// get point indices of corner vertices
	for ( i=0; i<8; ++i )
		corner[i] = grid_.point_idx( _cidx, i );

	// determine cube type
	for ( i=0; i<8; ++i )
		if ( grid_.scalar_distance( corner[i] ) > iso_ )
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


