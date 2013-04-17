/*:
 * OpenMeshX.h
 *
 *  Created on: Jul 16, 2012
 *      Author: hmirap
 */

#ifndef OPENMESHX_H_
#define OPENMESHX_H_

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <vector>
#include <utility>
#include <tuple>

struct MyTraits : public OpenMesh::DefaultTraits
{
  VertexAttributes(OpenMesh::Attributes::Status);
  FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color | OpenMesh::Attributes::Normal);
  EdgeAttributes(OpenMesh::Attributes::Status);
};


typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> OpenMeshExtended;

/*!
 * \struct OpenMeshXTraits
 * \brief traits of OpenMeshX
 * 
 * \ingroup OpenMeshX
 */
//template<typename TTraits>
//class mesh_traits< OpenMesh::PolyMesh_ArrayKernelT<TTraits> /*OpenMeshExtended*/>

template<>
class mesh_traits<OpenMeshExtended>

{
public:
	//typedef OpenMesh::PolyMesh_ArrayKernelT<TTraits> OpenMeshExtended;

	typedef OpenMesh::DefaultTraits original_traits;

	typedef typename OpenMeshExtended::Point Point;

	typedef typename OpenMeshExtended::Normal normal;


	typedef typename OpenMeshExtended::HalfedgeHandle h_edge_descriptor;//consider moving to advanced traits
	typedef typename OpenMeshExtended::HalfedgeIter h_edge_iterator;
	typedef typename std::vector<OpenMeshExtended::Halfedge>::size_type h_edges_size_type;

	typedef typename OpenMeshExtended::VertexHandle vertex_descriptor;
	typedef typename OpenMeshExtended::VertexIter vertex_iterator;
	typedef typename std::vector<OpenMeshExtended::Vertex>::size_type vertices_size_type;

	typedef typename OpenMeshExtended::EdgeHandle edge_descriptor;
	typedef typename OpenMeshExtended::EdgeIter edge_iterator;
	typedef typename std::vector<OpenMeshExtended::Edge>::size_type edges_size_type;

	typedef typename OpenMeshExtended::FaceHandle face_descriptor;
	typedef typename OpenMeshExtended::FaceIter face_iterator;
	typedef typename std::vector<OpenMeshExtended::Face>::size_type faces_size_type;

	class my_fv_iterator;
	class my_vv_iterator;
	class my_ve_iterator; 

	typedef my_fv_iterator fv_iterator;
	typedef my_vv_iterator vv_iterator;
	typedef my_ve_iterator ve_iterator;

    static std::pair<fv_iterator, fv_iterator>
    get_surrounding_vertices(const OpenMeshExtended& m_, face_descriptor fd);


//=================CONCEPTS======================


static bool remove_vertex(
				  vertex_descriptor v,
		  	  	  OpenMeshExtended &m);

static bool create_face(
				  vertex_descriptor a,
				  vertex_descriptor b,
				  vertex_descriptor c,
		  	  	  OpenMeshExtended& m);

static bool remove_face(
				  typename mesh_traits<OpenMeshExtended>::face_descriptor f,
		  	  	  OpenMeshExtended& m);

static std::pair<vertex_iterator,
	  	  	vertex_iterator>
get_all_vertices(const OpenMeshExtended& m_);

static std::pair<edge_iterator,
	  	  	edge_iterator>
get_all_edges(const OpenMeshExtended& m_);


static std::pair<face_iterator,
	  	  	face_iterator>
get_all_faces(const OpenMeshExtended& m_);

//=========== VERTEX ADJACENCY CONCEPT ===========

static bool is_isolated(const OpenMeshExtended& m_,
		vertex_descriptor v);

static std::pair<vv_iterator,
	  	  vv_iterator>
get_adjacent_vertices(
		const OpenMeshExtended& m_,
		  vertex_descriptor v);

static std::pair<ve_iterator, ve_iterator>
get_adjacent_edges(
		const OpenMeshExtended& m_,
		vertex_descriptor v);


 
};

template<>
class advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>> : public mesh_traits<OpenMeshExtended>
{
public:
	static mesh_traits<OpenMeshExtended>::normal 
	get_face_normal(
		const OpenMeshExtended& m_,
		face_descriptor f);

	template<
		typename Traits ,
		typename B = typename std::enable_if
			< Traits::FaceAttributes & OpenMesh::Attributes::Normal>::type
		>
	static bool 
	flip_face_normal_t(
		OpenMeshExtended& m_,
		face_descriptor& f)
        {

                typedef OpenMeshExtended Mesh;
                Mesh& m = const_cast<Mesh&>(m_);

                auto n = m.normal(f);
                m.set_normal(f, -n);
                return true;
	}


	static bool
	flip_face_normal(
		OpenMeshExtended& m_,
		face_descriptor& f)	{ return flip_face_normal_t<MyTraits>(m_, f); }

	static bool
	triangulate_face(
		OpenMeshExtended& m_,
		face_descriptor& f
	)
	{
                typedef OpenMeshExtended Mesh;
                Mesh& m = const_cast<Mesh&>(m_);

		m.triangulate(f);
	}

	static Point
	calc_middle_point(
		Point p1,
		Point p2,
		float coeff = 0.5f
	)
	{
		auto diff = p2;
		diff -= p1;
		diff *= coeff;
		diff += p2;
		return diff;
	}

	static vertex_descriptor
	calc_middle_vertex(
		OpenMeshExtended& m_,
		vertex_descriptor v1,
		vertex_descriptor v2,
		float coeff = 0.5f)
	{
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		auto diff = m.point( v2 );
		diff -= m.point( v1 );
		diff *= coeff;
		diff += m.point( v1 );


		return m.new_vertex( diff );
	}

	static bool
	split_edge(
		OpenMeshExtended& m_,
		edge_descriptor e)
	{
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		h_edge_descriptor heh     = m.halfedge_handle(e, 0);
		h_edge_descriptor opp_heh = m.halfedge_handle(e, 1);

		h_edge_descriptor new_heh, opp_new_heh, t_heh;
		vertex_descriptor   vh;
		vertex_descriptor   vh1( m.to_vertex_handle(heh));

		// new vertex
		vh                = calc_middle_vertex
			(m,
			m.to_vertex_handle(heh),
			m.to_vertex_handle(opp_heh));


		// Re-link mesh entities
		if (m.is_boundary(e))
		{
			for (t_heh = heh;
			m.next_halfedge_handle(t_heh) != opp_heh;
			t_heh = m.opposite_halfedge_handle(m.next_halfedge_handle(t_heh)))
			{}
		}
		else
		{
			for (t_heh = m.next_halfedge_handle(opp_heh);
			m.next_halfedge_handle(t_heh) != opp_heh;
			t_heh = m.next_halfedge_handle(t_heh) )
			{}
		}

		new_heh     = m.new_edge(vh, vh1);
		opp_new_heh = m.opposite_halfedge_handle(new_heh);
		m.set_vertex_handle( heh, vh );

		m.set_next_halfedge_handle(t_heh, opp_new_heh);
		m.set_next_halfedge_handle(new_heh, m.next_halfedge_handle(heh));
		m.set_next_halfedge_handle(heh, new_heh);
		m.set_next_halfedge_handle(opp_new_heh, opp_heh);


		if (m.face_handle(opp_heh).is_valid())
  		{
			m.set_face_handle(
				opp_new_heh,
				m.face_handle(opp_heh));
			m.set_halfedge_handle(
				m.face_handle(opp_new_heh), 
				opp_new_heh);
		}

		if( m.face_handle(heh).is_valid())
		{
			m.set_face_handle( new_heh, m.face_handle(heh) );
			m.set_halfedge_handle( m.face_handle(heh), heh );
		}

		m.set_halfedge_handle( vh, new_heh);
		m.set_halfedge_handle( vh1, opp_new_heh );

// Never forget this, when playing with the topology
		m.adjust_outgoing_halfedge( vh );
		m.adjust_outgoing_halfedge( vh1 );

		return true;
	}

	static bool
	truncate(
		OpenMeshExtended& m_,
		vertex_descriptor v,
		float coeff = 0.5f)
	{
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		typedef std::tuple<
			vertex_descriptor,	//vertex V
			h_edge_descriptor,	//HE1 incoming halfedge to V
			h_edge_descriptor,	//HE2 next halfedge directed by HE1
			h_edge_descriptor>	//HE3 opposite halfedge to HE1 pointing to V
				Vertex_Halfedge_Relation;

		typedef std::vector<Vertex_Halfedge_Relation> Relational_List;
		Relational_List r_l;


		auto he_start = m.halfedge_handle(v);

		auto he_1 = m.halfedge_handle(v);
		auto prev_he_1 = get_previous_halfedge(m, he_1);
		auto opp_he_1 = m.opposite_halfedge_handle(he_1);
		auto v_last = v;
		auto v_new = v;

		do
		{
			if (he_1 != he_start)
			v_new = calc_middle_vertex(
					m,
					v,
					m.to_vertex_handle(he_1),
					coeff);

			prev_he_1 = opp_he_1;
			he_1 = m.next_halfedge_handle(opp_he_1);
			opp_he_1 = m.opposite_halfedge_handle(he_1);

			r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));

		}
		while (he_1 != he_start);

		auto fin_center = calc_middle_point(
			m.point(v),
			m.point(m.to_vertex_handle(he_1)),
			-1.f * (1.f - coeff));

		m.set_point(v, fin_center);
		m.adjust_outgoing_halfedge(v);

		auto f_a = m.new_face();

		std::vector<h_edge_descriptor> inside_face_he;

		for (int i=0;i<r_l.size();i++)
		{
			auto vx_rel = r_l[i];
			auto vx_rel_b = r_l[(i+1) % r_l.size()];

			auto v_a = std::get<0>(vx_rel);
			auto v_b = std::get<0>(vx_rel_b);
			auto prev_he_a = std::get<1>(vx_rel);
			auto he_a = std::get<2>(vx_rel);
			auto opp_he_a = std::get<3>(vx_rel);

			auto he_new = m.new_edge(v_a, v_b);
			m.set_face_handle(m.opposite_halfedge_handle(he_new), f_a);
			m.set_halfedge_handle(f_a, m.opposite_halfedge_handle(he_new));

			inside_face_he.push_back(m.opposite_halfedge_handle(he_new));

			m.set_face_handle(he_new, m.face_handle(he_1));
			m.set_next_halfedge_handle(prev_he_a, he_new);
			m.set_next_halfedge_handle(he_new, he_a);
			m.set_vertex_handle(opp_he_a, v_b);
		
			m.set_halfedge_handle(v_a, he_new);
			m.adjust_outgoing_halfedge(v_a);
		}

		for (int i=0; i<inside_face_he.size(); i++)
		{
			auto hx_a = inside_face_he[ i ];
			auto hx_b = inside_face_he[ (i+1) % inside_face_he.size()];

			m.set_next_halfedge_handle(hx_b, hx_a);
		}

		return true;
	}	


	static bool
	bevel(
		OpenMeshExtended& m_,
		edge_descriptor e,
		float coeff = 0.5f)
	{
		
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		typedef std::tuple<
			vertex_descriptor,	//vertex V
			h_edge_descriptor,	//HE1 incoming halfedge to V
			h_edge_descriptor,	//HE2 next halfedge directed by HE1
			h_edge_descriptor>	//HE3 opposite halfedge to HE1 pointing to V
				Vertex_Halfedge_Relation;

		typedef std::vector<Vertex_Halfedge_Relation> Relational_List;
		Relational_List r_l;

		auto he_a = m.halfedge_handle(e, 0);
		auto he_b = m.halfedge_handle(e, 1);
		auto v_a = m.to_vertex_handle(he_a);
		auto v_b = m.to_vertex_handle(he_b);

		auto prev_he_1 = he_a;
		auto he_1 = m.next_halfedge_handle(prev_he_1);
		auto opp_he_1 = m.opposite_halfedge_handle(he_1);
		auto v_new = v_a;

		r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));

		do
		{
			prev_he_1 = opp_he_1;
			he_1 = m.next_halfedge_handle(opp_he_1);
			opp_he_1 = m.opposite_halfedge_handle(he_1);
			v_new = calc_middle_vertex(
					m,
					v_a,
					m.to_vertex_handle(he_1),
					coeff);

			r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));

		}
		while (m.next_halfedge_handle(opp_he_1) != he_b);

		prev_he_1 = he_b;
		he_1 = m.next_halfedge_handle(prev_he_1);
		opp_he_1 = m.opposite_halfedge_handle(he_1);
		v_new = calc_middle_vertex(
			m,
			v_b,
			m.to_vertex_handle(he_1),
			coeff);

		r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));

		do
		{	
			prev_he_1 = opp_he_1;
			he_1 = m.next_halfedge_handle(opp_he_1);
			opp_he_1 = m.opposite_halfedge_handle(he_1);

			if (m.next_halfedge_handle(opp_he_1) != he_a)
				v_new = calc_middle_vertex(
					m,
					v_b,
					m.to_vertex_handle(he_1),
					coeff);
			else
				v_new = v_b;


			r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));
		}
		while (m.next_halfedge_handle(opp_he_1) != he_a);

		auto f_a = m.new_face();

		std::vector<h_edge_descriptor> inside_face_he;

		for (int i=0; i<r_l.size() - 1; i++)
		{
			auto vx_rel = r_l[i];
			auto vx_rel_b = r_l[(i+1)];

			auto vx_a = std::get<0>(vx_rel);
			auto vx_b = std::get<0>(vx_rel_b);
			auto prev_hex_a = std::get<1>(vx_rel);
			auto prev_hex_b = std::get<1>(vx_rel_b);
			auto hex_a = std::get<2>(vx_rel);
			auto hex_b = std::get<2>(vx_rel_b);
			auto opp_hex_a = std::get<3>(vx_rel);
			auto opp_hex_b = std::get<3>(vx_rel_b);

			auto he_new = m.new_edge(vx_a, vx_b);

			m.set_face_handle(m.opposite_halfedge_handle(he_new), f_a);
			m.set_halfedge_handle(f_a, m.opposite_halfedge_handle(he_new));

			inside_face_he.push_back(m.opposite_halfedge_handle(he_new));

			m.set_face_handle(he_new, m.face_handle(opp_hex_a));
			m.set_next_halfedge_handle(opp_hex_a, he_new);
			m.set_next_halfedge_handle(he_new, hex_b);
			m.set_vertex_handle(opp_hex_b, vx_b);
		
			m.set_halfedge_handle(vx_a, he_new);
			m.adjust_outgoing_halfedge(vx_a);
		}

		m.set_face_handle(he_a, f_a);
		inside_face_he.push_back(he_b);

		for (int i=0; i<inside_face_he.size(); i++)
		{
			auto hx_a = inside_face_he[ i ];
			auto hx_b = inside_face_he[ (i+1) % inside_face_he.size()];
			m.set_next_halfedge_handle(hx_b, hx_a);
		}

		auto v_2a = m.to_vertex_handle(m.next_halfedge_handle(he_a));
		auto v_2b = m.from_vertex_handle(get_previous_halfedge(m, he_a));

		auto mid_p_a = calc_middle_point(m.point(v_a), m.point(v_2a), -1.f * (1.f - coeff));

		auto mid_p_b = calc_middle_point(m.point(v_b), m.point(v_2b), -1.f * (1.f - coeff));

		m.set_point(v_a, mid_p_a);
		m.adjust_outgoing_halfedge(v_a);

		m.set_point(v_b, mid_p_b);
		m.adjust_outgoing_halfedge(v_b);

		return true;
	}

	static
	h_edge_descriptor
	get_previous_halfedge(
		OpenMeshExtended& m,
		h_edge_descriptor heh)
	{
		auto t_heh = heh;

			for (t_heh = m.next_halfedge_handle(t_heh);
                        m.next_halfedge_handle(t_heh) != heh;
                        t_heh = m.next_halfedge_handle(t_heh) )
                        {}

		return t_heh;

	}

};



typedef mesh_traits<OpenMeshExtended> OpenMeshXTraits;

#include "OpenMeshX.tcc"
		  
#endif /* OPENMESHX_H_ */
