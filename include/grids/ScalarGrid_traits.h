#include "ScalarGridT.hh"

template < class F, template <class F > class GridT>
class ScalarGrid_traits
{
public:
	//typedef typename GridT<F>	Grid;
	typedef typename GridT<F>::PointIdx		Point_descriptor;
	typedef typename GridT<F>::CubeIdx		Cube_descriptor;
	typedef typename GridT<F>::CubeIterator	Cube_iterator;
	typedef typename GridT<F>::Vec3		Coordinates_descriptor;

	static
	bool is_inside( GridT<F> g, Point_descriptor p )	{return g.is_inside(p);}
	
	static
	int get_cube_corner(GridT<F> g, Cube_descriptor c, int i)	{return g.point_idx( c, i );}
	
	static
	Coordinates_descriptor get_coords(GridT<F> g, Point_descriptor p)	{return g.point( p );}
};