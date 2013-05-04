#include <geometry/ray_face_intersection.hpp>

template <
	typename TGrid,
	typename TMesh,
	typename TGrid_Traits,// = grid_traits<TGrid>, 
	typename TMesh_Traits = mesh_traits<TMesh>
	>
class Voxelize
{
public:
	typedef TGrid						Grid;

	typedef typename TGrid_Traits::Point_descriptor	Point_descriptor;
	typedef typename TGrid_Traits::Coordinates_descriptor Coordinates_descriptor;
	typedef typename TGrid_Traits::Point_scalar_type	Point_scalar_type;
	typedef typename TGrid_Traits::Point_properties	Point_properties;
	typedef typename TGrid_Traits::Cube_descriptor		Cube_descriptor;
	typedef typename TGrid_Traits::Cube_iterator		Cube_iterator;

	typedef typename TMesh_Traits::Point			Mesh_Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;

	
	Grid&	grid_;
	TMesh&            mesh_;
	
	~Voxelize(){}
	Voxelize( Grid& _grid, TMesh& _mesh)
		: grid_( _grid ),
		mesh_( _mesh )
	{
		for ( auto cube : grid_ )
			process_cube( cube );
		
		floodfill();
	}
	
	void process_cube(Cube_descriptor cube)
	{
		Point_descriptor	corner[8];
		Vertex_descriptor	samples[12];
		unsigned int		cubetype( 0 );
		unsigned int		i;

		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, cube, i );
		
		crosses_face( corner[0], corner[1] ,(Face_descriptor)0 );
		crosses_face( corner[1], corner[2] ,(Face_descriptor)0 );
		crosses_face( corner[3], corner[2] ,(Face_descriptor)0 );
		crosses_face( corner[0], corner[3] ,(Face_descriptor)0 );
		crosses_face( corner[4], corner[5] ,(Face_descriptor)0 );
		crosses_face( corner[5], corner[6] ,(Face_descriptor)0 );
		crosses_face( corner[7], corner[6] ,(Face_descriptor)0 );
		crosses_face( corner[4], corner[7] ,(Face_descriptor)0 );
		crosses_face( corner[0], corner[4] ,(Face_descriptor)0 );
		crosses_face( corner[1], corner[5] ,(Face_descriptor)0 );
		crosses_face( corner[2], corner[6] ,(Face_descriptor)0 );
		crosses_face( corner[3], corner[7] ,(Face_descriptor)0 );
				
		crosses_face( corner[0], corner[1] ,(Face_descriptor)1 );
		crosses_face( corner[1], corner[2] ,(Face_descriptor)1 );
		crosses_face( corner[3], corner[2] ,(Face_descriptor)1 );
		crosses_face( corner[0], corner[3] ,(Face_descriptor)1 );
		crosses_face( corner[4], corner[5] ,(Face_descriptor)1 );
		crosses_face( corner[5], corner[6] ,(Face_descriptor)1 );
		crosses_face( corner[7], corner[6] ,(Face_descriptor)1 );
		crosses_face( corner[4], corner[7] ,(Face_descriptor)1 );
		crosses_face( corner[0], corner[4] ,(Face_descriptor)1 );
		crosses_face( corner[1], corner[5] ,(Face_descriptor)1 );
		crosses_face( corner[2], corner[6] ,(Face_descriptor)1 );
		crosses_face( corner[3], corner[7] ,(Face_descriptor)1 );
				
		crosses_face( corner[0], corner[1] ,(Face_descriptor)2 );
		crosses_face( corner[1], corner[2] ,(Face_descriptor)2 );
		crosses_face( corner[3], corner[2] ,(Face_descriptor)2 );
		crosses_face( corner[0], corner[3] ,(Face_descriptor)2 );
		crosses_face( corner[4], corner[5] ,(Face_descriptor)2 );
		crosses_face( corner[5], corner[6] ,(Face_descriptor)2 );
		crosses_face( corner[7], corner[6] ,(Face_descriptor)2 );
		crosses_face( corner[4], corner[7] ,(Face_descriptor)2 );
		crosses_face( corner[0], corner[4] ,(Face_descriptor)2 );
		crosses_face( corner[1], corner[5] ,(Face_descriptor)2 );
		crosses_face( corner[2], corner[6] ,(Face_descriptor)2 );
		crosses_face( corner[3], corner[7] ,(Face_descriptor)2 );
						
		crosses_face( corner[0], corner[1] ,(Face_descriptor)3 );
		crosses_face( corner[1], corner[2] ,(Face_descriptor)3 );
		crosses_face( corner[3], corner[2] ,(Face_descriptor)3 );
		crosses_face( corner[0], corner[3] ,(Face_descriptor)3 );
		crosses_face( corner[4], corner[5] ,(Face_descriptor)3 );
		crosses_face( corner[5], corner[6] ,(Face_descriptor)3 );
		crosses_face( corner[7], corner[6] ,(Face_descriptor)3 );
		crosses_face( corner[4], corner[7] ,(Face_descriptor)3 );
		crosses_face( corner[0], corner[4] ,(Face_descriptor)3 );
		crosses_face( corner[1], corner[5] ,(Face_descriptor)3 );
		crosses_face( corner[2], corner[6] ,(Face_descriptor)3 );
		crosses_face( corner[3], corner[7] ,(Face_descriptor)3 );
	}
	
	bool crosses_face(Point_descriptor corner0, Point_descriptor corner1, Face_descriptor f)
	{
		auto c0 = TGrid_Traits::get_coords(grid_, corner0);
		auto c1 = TGrid_Traits::get_coords(grid_, corner1);
		auto fv_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		Vertex_descriptor v[3];
		
		auto fv = fv_pair.first;
		for(int i = 0; i < 3; i++ )
		{
			v[i] = *fv;
			++fv;
		}
		
		auto a = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[0]);
		auto b = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[1]);
		auto c = (Coordinates_descriptor)TMesh_Traits::get_coordinates(mesh_, v[2]);
		
		auto c_direction = c1 - c0;
		
		Coordinates_descriptor intersection;
		float distance = 0.f;
		
		auto result = line_face_intersection(
		c0,
		c_direction,
		a,
		b,
		c,
		intersection,
		distance);
		
		if ((c0 - c1).length() < distance)
			result = false;
		
// 		std::cout << "intersection: " << result << " & " << (c0 - c1).length() << " & " << distance << " & " << intersection << std::endl;
		
		if (result)
		{
			auto normal = TMesh_Traits::get_face_normal(mesh_, f);
			if (dot(normal, c_direction) < 0)
			{
				TGrid_Traits::set_scalar_value( grid_, corner0, 0.5);
				if (TGrid_Traits::scalar_value( grid_, corner1) == -0.5)
				{
					TGrid_Traits::set_scalar_value( grid_, corner1, -0.5);
				}
				else
				{
					TGrid_Traits::set_scalar_value( grid_, corner1, -0.5);
				}
			}
			else
			{
				TGrid_Traits::set_scalar_value( grid_, corner1, 0.5);
				if (TGrid_Traits::scalar_value( grid_, corner0) == -0.5)
				{
					TGrid_Traits::set_scalar_value( grid_, corner0, -0.5);
				}
				else
				{
					TGrid_Traits::set_scalar_value( grid_, corner0, -0.5);
				}
			}
		}
		return result;
	}
	
	void floodfill()
	{
		auto x = grid_.x_resolution();
		int i = 0;
		bool inside = false;
		bool fill = false;
		bool previous_boundary = false;
		std::vector<int> to_fill;
	
		for (auto cube: grid_)
		{
			if (i == x)
			{
				if (inside)
				{
					
				}
				else
				{
					for (auto c1 : to_fill)
					{
						fill_cube(c1);
					}
				}
				to_fill.clear();
				i = 0;
				inside = false;
				std::cout << "\n" << std::endl;
			}
			
			if (boundary_cube(cube))
			{
				if (inside && !previous_boundary)
				{
					inside = false;
				}
				else if (!previous_boundary)
				{
					inside = true;
				}
				
				previous_boundary = true;
			}
			else
			{
				previous_boundary = false;
			}
			/*
			if (inside && !boundary_cube(cube))
			{
				fill_cube(cube);
			}*/
			
			std::cout << "bound: " << boundary_cube(cube) << " io: " << (inside ? "in" : "out") << std::endl;
			if (!boundary_cube(cube) && inside)
			{
				to_fill.push_back(cube);
			}
			i++;
		}
	}
	
	void fill_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		TGrid_Traits::set_scalar_value( grid_, corner0, 0.5);
		TGrid_Traits::set_scalar_value( grid_, corner1, 0.5);
		
		return;
		for (int i=0; i<8; ++i )
		{
			auto corner = TGrid_Traits::get_cube_corner(grid_, _cidx, i );
			TGrid_Traits::set_scalar_value( grid_, corner, -0.5);
		}
	}
	
	bool boundary_cube(unsigned int _cidx)
	{
		Point_descriptor	corner[8];
		unsigned int		cubetype( 0 );
		int i;

		bool opp_bound = false;
		
		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, _cidx, i );

		return TGrid_Traits::is_inside( grid_, corner[0] ) ^ TGrid_Traits::is_inside( grid_, corner[1] );
		
		for (i=0; i<8; ++i )
		{
			if (TGrid_Traits::scalar_value( grid_, corner[i]) == -0.5)
			{
				opp_bound = true;
			}
			if ( !TGrid_Traits::is_inside( grid_, corner[i] ))
			{
				cubetype |= ( 1<<i );
			}
		}

		if ( cubetype == 0 || cubetype == 255 )
			return false;
		else
			return true && opp_bound;
	}
};