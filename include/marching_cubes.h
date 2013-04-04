#include <grids/Grid.hh>

template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
class MarchingCubes
{

	extern int edgeTable[256];
	extern int triTable[256][2][17];

	typedef Grid::PointIdx      PointIdx;
	typedef Grid::CubeIdx       CubeIdx;
	typedef Grid::CubeIterator  CubeIterator;

	typedef typename TMesh::Point         Point;
	typedef typename TMesh::VertexHandle  VertexHandle;
	typedef typename TMesh::FaceHandle    FaceHandle;

	template <class TMesh, class TMesh_Traits>
	int (IsoEx::Grid &g, TMesh& m)
	{

	}

	OpenMesh::VertexHandle
	add_vertex( PointIdx _p0, PointIdx _p1, IsoEx::ScalarGridT<double> grid_, OpenMeshExtended mesh_, float iso_, Edge2VertexMapT<PointIdx, VertexHandle> edge2vertex_)
	{
		// find vertex if it has been computed already
		VertexHandle   vh = edge2vertex_.find( _p0, _p1 );
		if ( vh.is_valid() )  return vh;
		// generate new vertex
		const OpenMeshExtended::Point&  p0( grid_.point( _p0 ) );
		const OpenMeshExtended::Point&  p1( grid_.point( _p1 ) );
		
	
		float s0 = grid_.scalar_distance( _p0 );
		float s1 = grid_.scalar_distance( _p1 );

		if ( fabs( s1-s0 ) > 0.00000001 )
		{
			float t  = ( iso_-s0 ) / ( s1-s0 );
			//vh = mesh_.add_vertex( OpenMesh::Vec3f(p0 + ( p1-p0 )*t) );
			vh = mesh_.add_vertex( OpenMeshExtended::Point((p0 + ( p1-p0 )*t)));

			return vh;
		}
		else
		{
			std::cout << s0 << ", " <<s1 <<"!!!!!!!!!!!!!!!!!" <<( s1 - s0 ) << std::endl;
			return vh;
		}
	}


	template <class Mesh>
	void
	MarchingCubesT<Mesh>::
	process_cube( CubeIdx _cidx )
	{

	PointIdx           corner[8];
	VertexHandle       samples[12];
	unsigned char      cubetype( 0 );
	unsigned int       i;


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
      mesh_.add_face( samples[triTable[cubetype][0][i  ]],
                      samples[triTable[cubetype][0][i+1]],
                      samples[triTable[cubetype][0][i+2]] );
	}
};


template <class PointIdx, class VertexHandle>
class Edge2VertexMapT
{
public:
   
  /// Constructor
  Edge2VertexMapT() {}


  /// clear the map
  void clear() { map_.clear(); }

  /// Store vertex in map
  void insert(PointIdx _p0, PointIdx _p1, VertexHandle _vhnd)
  {
    map_[EdgeKey(_p0, _p1)] = _vhnd;
  }

  /// Get vertex handle from map. Returns invalid handle if not found.
  VertexHandle find(PointIdx _p0, PointIdx _p1) const 
  {
    MyMapIterator it = map_.find(EdgeKey(_p0, _p1));
    if (it != map_.end())  return it->second;
    else return VertexHandle();
  }


private:

  class EdgeKey 
  {
  public:

    EdgeKey(PointIdx _p0, PointIdx _p1) {
      if (_p0 < _p1)  { p0_ = _p0;  p1_ = _p1; }
      else            { p0_ = _p1;  p1_ = _p0; }
    }

    bool operator<(const EdgeKey& _rhs) const 
    {
      if (p0_ != _rhs.p0_)
	return (p0_ < _rhs.p0_);
      else
	return (p1_ < _rhs.p1_);
    }

  private:
    PointIdx p0_, p1_;
  };


  typedef std::map<EdgeKey, VertexHandle>  MyMap;
  typedef typename MyMap::const_iterator   MyMapIterator;

  MyMap  map_;
};
