#include <tuple>

#include "ScalarGridT.hh"

template < class F, template <class F > class GridT>
class ScalarGrid_traits
{
public:
	//typedef typename GridT<F>	Grid;
	typedef typename GridT<F>::PointIdx		Point_descriptor;
	typedef F					Point_scalar_type;
	typedef std::tuple<float, int>			Point_properties;
	typedef typename GridT<F>::CubeIdx		Cube_descriptor;
	typedef typename GridT<F>::CubeIterator	Cube_iterator;
	typedef typename GridT<F>::Vec3		Coordinates_descriptor;

	static inline
	Point_scalar_type
	scalar_value( GridT<F> g, Point_descriptor p )	{return g.scalar_distance(p);}

	static inline
	void
	set_scalar_value( GridT<F>& g, Point_descriptor p, Point_scalar_type val )	{ g.set_scalar_distance(p, val);}

	static inline
	bool is_inside( GridT<F> g, Point_descriptor p )	{return g.scalar_distance(p)==0.5 /*|| g.scalar_distance(p)==0*/ ; /*g.is_inside(p);*/}
	
	static inline
	int get_cube_corner(GridT<F> g, Cube_descriptor c, int i)	{return g.point_idx( c, i );}
	
	static inline
	Coordinates_descriptor get_coords(GridT<F> g, Point_descriptor p)	{return g.point( p );}

	static inline
	Point_properties
	get_point_properties(GridT<F> g, Point_descriptor p )
	{
		return std::make_tuple( g.scalar_distance(p), 555);
	}

};
