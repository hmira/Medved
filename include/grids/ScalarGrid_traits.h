#include <tuple>
#include <boost/concept_check.hpp>

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
	typedef typename GridT<F>::Vec3		Coordinates_type;
	typedef float					Vector_unit_type;

	static inline
	Point_scalar_type
	scalar_value( GridT<F> g, Point_descriptor p )	{return g.scalar_distance(p);}

	static inline
	void
	set_scalar_value( GridT<F>& g, Point_descriptor p, Point_scalar_type val )	{ g.set_scalar_distance(p, val);}

	static inline
	void
	fill_cube_by_value( GridT<F>& g, Cube_descriptor c, Point_scalar_type val)
	{
		for (int i=0; i<8; ++i )
		{
			auto corner = g.point_idx( c, i );
			g.set_scalar_distance(corner, val);
		}
	}

	static inline
	bool
	is_cube_inside( GridT<F>& g, Cube_descriptor c)
	{
		for (int i=0; i<8; ++i )
		{
			auto corner = g.point_idx( c, i );
			if (g.scalar_distance(corner) <= 0)
				return false;
		}
		return true;
	}
	
	static inline
	bool is_inside( GridT<F> g, Point_descriptor p )	{return g.scalar_distance(p)==0.5 /*|| g.scalar_distance(p)==0*//* g.is_inside(p)*/;}
	
	static inline
	int get_cube_corner(GridT<F> g, Cube_descriptor c, int i)	{return g.point_idx( c, i );}
	
	static inline
	Coordinates_type get_coords(GridT<F> g, Point_descriptor p)	{return g.point( p );}

	static inline
	Point_properties
	get_point_properties(GridT<F> g, Point_descriptor p )
	{
		return std::make_tuple( g.scalar_distance(p), 555);
	}
	
	static inline
	std::pair<int, int>
	get_bounds(const GridT<F> &g, int num_of_ths, int ord)
	{
		auto x = g.x_resolution();
		auto y = g.y_resolution();
		auto z = g.z_resolution();
		
		auto len = x*y*z;
		auto low_bound = ord * ( x / num_of_ths ) * y * z;
		auto up_bound = (ord + 1) * ( x / num_of_ths - 2 ) * y * z;
		
		return std::make_pair(low_bound, up_bound);
	}
	
	static inline
	float
	get_x_element_size(const GridT<F> &g)
	{
		return g.dx().norm();
// 		return g.x_axis().norm() / g.x_resolution();
	}

	static inline
	float
	get_y_element_size(const GridT<F> &g)
	{
		return g.dy().norm();
// 		return g.y_axis().norm() / g.y_resolution();
	}

	static inline
	float
	get_z_element_size(const GridT<F> &g)
	{
		return g.dz().norm();
// 		return g.z_axis().norm() / g.z_resolution();
	}
	
	static inline
	int
	get_x_resolution(const GridT<F> &g)
	{
		return g.x_resolution() - 1;
	}

	static inline
	int
	get_y_resolution(const GridT<F> &g)
	{
		return g.y_resolution() - 1;
	}

	static inline
	int
	get_z_resolution(const GridT<F> &g)
	{
		return g.z_resolution() - 1;
	}
	
	static inline
	Cube_descriptor
	get_cube_from_coords(const GridT<F> &g, const Coordinates_type& pt)
	{
		unsigned int X( g.x_resolution() -1 ), Y( g.y_resolution() -1 );

// 		Vec3 pl = const_cast<Vec3&>( pt );
		
		auto dxf = g.dx().norm();
		auto dyf = g.dy().norm();
		auto dzf = g.dz().norm();
		
		int x = ( int )( pt[0] / dxf );
		int y = ( int )( pt[1] / dyf );
		int z = ( int )( pt[2] / dzf );

		return x + y * X + z * X * Y;
	}

};
