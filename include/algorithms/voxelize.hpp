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
	typedef typename TGrid_Traits::Vector_unit_type	Grid_unit_type;

	typedef typename TMesh_Traits::Point			Mesh_Point;
	typedef typename TMesh_Traits::vertex_descriptor	Vertex_descriptor;
	typedef typename TMesh_Traits::face_descriptor		Face_descriptor;

	
	Grid&	grid_;
	Grid_unit_type elem_x_, elem_y_, elem_z_;
	TMesh&            mesh_;
	
	~Voxelize(){}
	Voxelize( Grid& _grid, TMesh& _mesh)
		: grid_( _grid ),
		mesh_( _mesh )
	{
		elem_x_ = TGrid_Traits::get_x_element_size(grid_);
		elem_y_ = TGrid_Traits::get_y_element_size(grid_);
		elem_z_ = TGrid_Traits::get_z_element_size(grid_);
		
		auto all_faces = TMesh_Traits::get_all_faces(mesh_);
		for (auto f = all_faces.first; f != all_faces.second; ++f)
		{
			Vertex_descriptor vertices[3];
			int i = 0;
			auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
			for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
			{
				vertices[i] = *fv_i;
				i++;
			}
			
			auto cp0 = TMesh_Traits::get_coordinates(mesh_, vertices[0]);
			auto cp1 = TMesh_Traits::get_coordinates(mesh_, vertices[1]);
			auto cp2 = TMesh_Traits::get_coordinates(mesh_, vertices[2]);

			if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
			if (cp1[1] > cp2[1]) std::swap(cp1, cp2);
			if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
			
			auto minY = cp0[1];
			auto maxY = cp2[1];
			
			if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
			if (cp1[0] > cp2[0]) std::swap(cp1, cp2);
			if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
			
			auto minX = cp0[0];
			auto maxX = cp2[0];
			
			if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
			if (cp1[2] > cp2[2]) std::swap(cp1, cp2);
			if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
			
			auto minZ = cp0[2];
			auto maxZ = cp2[2];
			
			/*try*/
			
			auto dx = maxX - minX;
			auto dy = maxY - minY;
			auto dz = maxZ - minZ;
			
			rasterize_triangle(*f);
			
			if (fabs(dz) < fabs(dx) && fabs(dz) < fabs(dy))
				fill_triangle_xy(*f);
			else if (fabs(dy) < fabs(dx) && fabs(dy) < fabs(dz))
				fill_triangle_xz(*f);
			else
				fill_triangle_yz(*f);
		}

// 		for ( auto cube : grid_ )
// 			process_cube( cube );
// 		
// 		floodfill();
	}
	
	void process_cube(Cube_descriptor cube)
	{
		Point_descriptor	corner[8];
		Vertex_descriptor	samples[12];
		unsigned int		cubetype( 0 );
		unsigned int		i;

		for ( i=0; i<8; ++i )
			corner[i] = TGrid_Traits::get_cube_corner(grid_, cube, i );

		
		
// 		crosses_face( corner[0], corner[1] ,(Face_descriptor)0 );
// 		crosses_face( corner[1], corner[2] ,(Face_descriptor)0 );
// 		crosses_face( corner[3], corner[2] ,(Face_descriptor)0 );
// 		crosses_face( corner[0], corner[3] ,(Face_descriptor)0 );
// 		crosses_face( corner[4], corner[5] ,(Face_descriptor)0 );
// 		crosses_face( corner[5], corner[6] ,(Face_descriptor)0 );
// 		crosses_face( corner[7], corner[6] ,(Face_descriptor)0 );
// 		crosses_face( corner[4], corner[7] ,(Face_descriptor)0 );
// 		crosses_face( corner[0], corner[4] ,(Face_descriptor)0 );
// 		crosses_face( corner[1], corner[5] ,(Face_descriptor)0 );
// 		crosses_face( corner[2], corner[6] ,(Face_descriptor)0 );
// 		crosses_face( corner[3], corner[7] ,(Face_descriptor)0 );
// 				
// 		crosses_face( corner[0], corner[1] ,(Face_descriptor)1 );
// 		crosses_face( corner[1], corner[2] ,(Face_descriptor)1 );
// 		crosses_face( corner[3], corner[2] ,(Face_descriptor)1 );
// 		crosses_face( corner[0], corner[3] ,(Face_descriptor)1 );
// 		crosses_face( corner[4], corner[5] ,(Face_descriptor)1 );
// 		crosses_face( corner[5], corner[6] ,(Face_descriptor)1 );
// 		crosses_face( corner[7], corner[6] ,(Face_descriptor)1 );
// 		crosses_face( corner[4], corner[7] ,(Face_descriptor)1 );
// 		crosses_face( corner[0], corner[4] ,(Face_descriptor)1 );
// 		crosses_face( corner[1], corner[5] ,(Face_descriptor)1 );
// 		crosses_face( corner[2], corner[6] ,(Face_descriptor)1 );
// 		crosses_face( corner[3], corner[7] ,(Face_descriptor)1 );
// 				
// 		crosses_face( corner[0], corner[1] ,(Face_descriptor)2 );
// 		crosses_face( corner[1], corner[2] ,(Face_descriptor)2 );
// 		crosses_face( corner[3], corner[2] ,(Face_descriptor)2 );
// 		crosses_face( corner[0], corner[3] ,(Face_descriptor)2 );
// 		crosses_face( corner[4], corner[5] ,(Face_descriptor)2 );
// 		crosses_face( corner[5], corner[6] ,(Face_descriptor)2 );
// 		crosses_face( corner[7], corner[6] ,(Face_descriptor)2 );
// 		crosses_face( corner[4], corner[7] ,(Face_descriptor)2 );
// 		crosses_face( corner[0], corner[4] ,(Face_descriptor)2 );
// 		crosses_face( corner[1], corner[5] ,(Face_descriptor)2 );
// 		crosses_face( corner[2], corner[6] ,(Face_descriptor)2 );
// 		crosses_face( corner[3], corner[7] ,(Face_descriptor)2 );
// 						
// 		crosses_face( corner[0], corner[1] ,(Face_descriptor)3 );
// 		crosses_face( corner[1], corner[2] ,(Face_descriptor)3 );
// 		crosses_face( corner[3], corner[2] ,(Face_descriptor)3 );
// 		crosses_face( corner[0], corner[3] ,(Face_descriptor)3 );
// 		crosses_face( corner[4], corner[5] ,(Face_descriptor)3 );
// 		crosses_face( corner[5], corner[6] ,(Face_descriptor)3 );
// 		crosses_face( corner[7], corner[6] ,(Face_descriptor)3 );
// 		crosses_face( corner[4], corner[7] ,(Face_descriptor)3 );
// 		crosses_face( corner[0], corner[4] ,(Face_descriptor)3 );
// 		crosses_face( corner[1], corner[5] ,(Face_descriptor)3 );
// 		crosses_face( corner[2], corner[6] ,(Face_descriptor)3 );
// 		crosses_face( corner[3], corner[7] ,(Face_descriptor)3 );
	}
	
	bool rasterize_triangle(Face_descriptor f)
	{
		auto f_edge_pair = TMesh_Traits::get_surrounding_edges(mesh_, f);
		for (auto fe_i = f_edge_pair.first; fe_i != f_edge_pair.second; ++fe_i)
		{
			auto v_pair = TMesh_Traits::get_edge_vertices(mesh_, *fe_i);
			auto v0 = v_pair.first;
			auto v1 = v_pair.second;
			
			std::cout << "pair: " << v0 << " : " << v1 << std::endl;
			rasterize_line(v0, v1);
			
		}
		
		
	}
	
	bool rasterize_line(Vertex_descriptor v0, Vertex_descriptor v1)
	{
		auto p0 = TMesh_Traits::get_coordinates(mesh_, v0);
		auto p1 = TMesh_Traits::get_coordinates(mesh_, v1);
		std::cout << p0[0] << p0[1] << p0[2] << std::endl;
		
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		if (fabs(dx) > fabs(dy) && fabs(dx) > fabs(dz))
		{
			rasterize_line_x(p0, p1);
		}
		
		else if (fabs(dy) > fabs(dx) && fabs(dy) > fabs(dz))
		{
			rasterize_line_y(p0, p1);
		}
		
		else if (fabs(dz) > fabs(dx) && fabs(dz) > fabs(dy))
		{
			rasterize_line_z(p0, p1);
		}
	}
	
	bool rasterize_line_x(Mesh_Point &p0, Mesh_Point &p1)
	{	
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0x < p1x) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dyx = dy / dx;
		auto dzx = dz / dx;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( ( x < p1x - elem_x_ ) || (x > p1x + elem_x_ ))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = grid_.nearest_point( v );
			TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
			
			x += ninja_constant * elem_x_;
			y += ninja_constant * dyx * elem_x_;
			z += ninja_constant * dzx * elem_x_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = grid_.nearest_point( v );
		TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
	}
	bool rasterize_line_y(Mesh_Point &p0, Mesh_Point &p1)
	{		
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0y < p1y) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dxy = dx / dy;
		auto dzy = dz / dy;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( (y < p1y - elem_y_) || (y > p1y + elem_y_ ))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = grid_.nearest_point( v );
			TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
			
			x += ninja_constant * dxy * elem_y_;
			y += ninja_constant * elem_y_;
			z += ninja_constant * dzy * elem_y_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = grid_.nearest_point( v );
		TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
	}
	bool rasterize_line_z(Mesh_Point &p0, Mesh_Point &p1)
	{
		auto p0x = p0[0];
		auto p0y = p0[1];
		auto p0z = p0[2];
		
		auto p1x = p1[0];
		auto p1y = p1[1];
		auto p1z = p1[2];
		
		int ninja_constant = (p0z < p1z) ? 1 : -1;
		
		auto dx = p0x - p1x;
		auto dy = p0y - p1y;
		auto dz = p0z - p1z;
		
		auto dxz = dx / dz;
		auto dyz = dy / dz;
		
		auto x = p0x;
		auto y = p0y;
		auto z = p0z;
		
		while ( (z < p1z - elem_z_) || (z > p1z + elem_z_))
		{
			Coordinates_descriptor v(x,y,z);
			auto pt = grid_.nearest_point( v );
			TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
			
			x += ninja_constant * dxz * elem_z_;
			y += ninja_constant * dyz * elem_z_;
			z += ninja_constant * elem_z_;
		}
		
		Coordinates_descriptor v(x,y,z);
		auto pt = grid_.nearest_point( v );
		TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
	}
	
	bool fill_triangle_xy(Face_descriptor f)
	{
		Vertex_descriptor vertices[3];
		int i = 0;
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			vertices[i] = *fv_i;
			i++;
		}
		
		auto cp0 = TMesh_Traits::get_coordinates(mesh_, vertices[0]);
		auto cp1 = TMesh_Traits::get_coordinates(mesh_, vertices[1]);
		auto cp2 = TMesh_Traits::get_coordinates(mesh_, vertices[2]);

		auto v0 = (cp1 - cp0);
		auto v1 = (cp2 - cp0);
		auto norm = cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * cp0[0] + b * cp0[1] + c * cp0[2]);
		
		if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
		if (cp1[1] > cp2[1]) std::swap(cp1, cp2);
		if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
		
		auto minY = cp0[1];
		auto maxY = cp2[1];
		
		if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
		if (cp1[0] > cp2[0]) std::swap(cp1, cp2);
		if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
		
		auto minX = cp0[0];
		auto maxX = cp2[0];
		
		/*try*/
		
		auto dx = maxX - minX;
		auto dy = maxY - minY;

		auto dyx = dy / dx;
		auto dxy = dx / dy;
		
		auto x = minX;
		auto y = minY;
		
		bool parallel_to_axis = (c < 0.00001 && c > -0.00001);
		
		if (parallel_to_axis) std::cout << "si v piči: " << norm << " " << d << std::endl;
		
// 		if (c == 0) return false;
		auto rcp_c = 1.f / c;
		auto z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
		
		if (dx > dy)
		{
			while (( y < maxY - dyx * elem_x_) || (y > maxY + dyx * elem_x_))
			{
				auto x = minX;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				
				while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					x += elem_x_;
					z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				}
				y += dyx * elem_x_;
			}
		}
		else
		{
			while (( x < maxX - dxy * elem_y_) || (x > maxX + dxy * elem_y_))
			{
				y = minY;
				z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				
				while (( y < maxY - elem_y_ ) || (y > maxY + elem_y_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					y += elem_y_;
					z = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + b * y) * rcp_c;
				}
				x += dxy * elem_y_;
			}
		}
	}
	
	bool fill_triangle_xz(Face_descriptor f)
	{
		Vertex_descriptor vertices[3];
		int i = 0;
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			vertices[i] = *fv_i;
			i++;
		}
		
		auto cp0 = TMesh_Traits::get_coordinates(mesh_, vertices[0]);
		auto cp1 = TMesh_Traits::get_coordinates(mesh_, vertices[1]);
		auto cp2 = TMesh_Traits::get_coordinates(mesh_, vertices[2]);

		auto v0 = (cp1 - cp0);
		auto v1 = (cp2 - cp0);
		auto norm = cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * cp0[0] + b * cp0[1] + c * cp0[2]);
		
		if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
		if (cp1[2] > cp2[2]) std::swap(cp1, cp2);
		if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
		
		auto minZ = cp0[2];
		auto maxZ = cp2[2];
		
		if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
		if (cp1[0] > cp2[0]) std::swap(cp1, cp2);
		if (cp0[0] > cp1[0]) std::swap(cp0, cp1);
		
		auto minX = cp0[0];
		auto maxX = cp2[0];
		
		/*try*/
		
		auto dx = maxX - minX;
		auto dz = maxZ - minZ;

		auto dzx = dz / dx;
		auto dxz = dx / dz;
		
		auto x = minX;
		auto z = minZ;
		
		bool parallel_to_axis = (b < 0.00001 && b > -0.00001);
		auto rcp_b = 1.f / b;
		auto y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
		
		if (dx > dz)
		{
			while (( z < maxZ - dzx * elem_x_) || (z > maxZ + dzx * elem_x_))
			{
				auto x = minX;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				
				while (( x < maxX - elem_x_ ) || (x > maxX + elem_x_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					x += elem_x_;
					y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				}
				z += dzx * elem_x_;
			}
		}
		else
		{
			while (( x < maxX - dxz * elem_z_) || (x > maxX + dxz * elem_z_))
			{
				auto z = minZ;
				y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				
				while (( z < maxZ - elem_z_ ) || (z > maxZ + elem_z_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					z += elem_z_;
					y = (parallel_to_axis) ? -1 * d : -1 * (d + a * x + c * z) * rcp_b;
				}
				x += dxz * elem_z_;
			}
		}
	}
	
	bool fill_triangle_yz(Face_descriptor f)
	{
		Vertex_descriptor vertices[3];
		int i = 0;
		auto f_vertex_pair = TMesh_Traits::get_surrounding_vertices(mesh_, f);
		for (auto fv_i = f_vertex_pair.first; fv_i != f_vertex_pair.second; ++fv_i)
		{
			vertices[i] = *fv_i;
			i++;
		}
		
		auto cp0 = TMesh_Traits::get_coordinates(mesh_, vertices[0]);
		auto cp1 = TMesh_Traits::get_coordinates(mesh_, vertices[1]);
		auto cp2 = TMesh_Traits::get_coordinates(mesh_, vertices[2]);

		auto v0 = (cp1 - cp0);
		auto v1 = (cp2 - cp0);
		auto norm = cross(v0, v1);
		
		/*parameters*/
		auto a = norm[0];
		auto b = norm[1];
		auto c = norm[2];
		auto d = -1 * (a * cp0[0] + b * cp0[1] + c * cp0[2]);
		
		if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
		if (cp1[1] > cp2[1]) std::swap(cp1, cp2);
		if (cp0[1] > cp1[1]) std::swap(cp0, cp1);
		
		auto minY = cp0[1];
		auto maxY = cp2[1];
		
		if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
		if (cp1[2] > cp2[2]) std::swap(cp1, cp2);
		if (cp0[2] > cp1[2]) std::swap(cp0, cp1);
		
		auto minZ = cp0[2];
		auto maxZ = cp2[2];
		
		/*try*/
		
		auto dz = maxZ - minZ;
		auto dy = maxY - minY;

		auto dyz = dy / dz;
		auto dzy = dz / dy;
		
		auto z = minZ;
		auto y = minY;
		bool parallel_to_axis = (a < 0.00001 && a > -0.00001);
		auto rcp_a = 1.f / a;
		auto x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
		
		if (dz > dy)
		{
			while (( y < maxY - dyz * elem_z_) || (y > maxY + dyz * elem_z_))
			{
				auto z = minZ;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				
				while (( z < maxZ - elem_z_ ) || ( z > maxZ + elem_z_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					z += elem_z_;
					x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				}
				y += dyz * elem_z_;
			}
		}
		else
		{
			while (( z < maxZ - dzy * elem_y_) || (z > maxZ + dzy * elem_y_))
			{
				auto y = minY;
				x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				
				while (( y < maxY - elem_y_ ) || ( y > maxY + elem_y_ ))
				{
					if ((x > 0) && (y > 0) && (z > 0))
					{
						Coordinates_descriptor v(x,y,z);
						auto pt = grid_.nearest_point( v );
						TGrid_Traits::set_scalar_value(grid_, pt, 0.5);
					}
					y += elem_y_;
					x = (parallel_to_axis) ? -1 * d : -1 * (d + c * z + b * y) * rcp_a;
				}
				z += dzy * elem_y_;
			}
		}
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
			if (dot(normal, c_direction) > 0)
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
			
			std::cout << "bound: " << boundary_cube(cube) << " io: " << (inside ? "in" : "out")
			<< " start: " << (starting_cube(cube) ? " y " : " n ") << " end: " << ( ending_cube(cube) ? " y " : " n " )  << std::endl;
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
	
	bool starting_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		auto a = TGrid_Traits::scalar_value( grid_, corner0);
		auto b = TGrid_Traits::scalar_value( grid_, corner1);
		
		return (a == -0.5) && (b == 0.5);
		
		return !TGrid_Traits::is_inside( grid_, corner0 ) && TGrid_Traits::is_inside( grid_, corner1 );
	}
	
	bool ending_cube(unsigned int _cidx)
	{
		auto corner0 = TGrid_Traits::get_cube_corner(grid_, _cidx, 0 );
		auto corner1 = TGrid_Traits::get_cube_corner(grid_, _cidx, 1 );
		auto a = TGrid_Traits::scalar_value( grid_, corner0);
		auto b = TGrid_Traits::scalar_value( grid_, corner1);
		
		return (a == 0.5) && (b == -0.5);
		
		return TGrid_Traits::is_inside( grid_, corner0 ) && !TGrid_Traits::is_inside( grid_, corner1 );
	}
	
	bool boundary_cube(unsigned int _cidx)
	{
		return starting_cube(_cidx) || ending_cube(_cidx);
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