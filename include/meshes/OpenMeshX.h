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

	static vertex_descriptor
	calc_middle_vertex(
		OpenMeshExtended& m_,
		vertex_descriptor v1,
		vertex_descriptor v2)
	{
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		auto middle_point = m.point( v1 );
		middle_point += m.point( v2 );
		middle_point *= 0.5;

		return m.new_vertex( middle_point );
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
		vertex_descriptor v)
	{
		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		auto he_1 = m.halfedge_handle(v);
		auto he_2 = m.opposite_halfedge_handle(he_1);
		auto he_3 = m.next_halfedge_handle(he_2);
		auto he_4 = m.opposite_halfedge_handle(he_3);
		auto he_5 = m.next_halfedge_handle(he_4);
		auto he_6 = m.opposite_halfedge_handle(he_5);

		if (m.next_halfedge_handle(he_6) != he_1)
			return false;
/*
		std::cout
		<< he_1 << ":1\n" 
		<< he_2 << ":2\n"
		<< he_3 << ":3\n"
		<< he_4 << ":4\n"
		<< he_5 << ":5\n"
		<< he_6 << ":6\n" 
		<< m.next_halfedge_handle(he_6) << ":6_opp\n" << std::endl;
*/

		auto v1 = calc_middle_vertex(
			m,
			v,
			m.to_vertex_handle(he_1));

		auto v2 = calc_middle_vertex(
			m,
			v,
			m.to_vertex_handle(he_3));

		m.set_vertex_handle(he_2, v1);
		m.set_vertex_handle(he_4, v2);

		auto new_he_1 = m.new_edge(v, v1);
		auto new_he_2 = m.new_edge(v1, v2);
		auto new_he_3 = m.new_edge(v2, v);

		m.set_halfedge_handle(v1, he_1);
		m.set_halfedge_handle(v2, he_3);

		m.set_next_halfedge_handle(he_6, new_he_1);
		m.set_next_halfedge_handle(new_he_1, he_1);
		m.set_next_halfedge_handle(he_2, new_he_2);
		m.set_next_halfedge_handle(new_he_2, he_3);
		m.set_next_halfedge_handle(he_4, new_he_3);
		m.set_next_halfedge_handle(new_he_3, he_6);


		m.set_face_handle( new_he_1, m.face_handle(he_1));
		m.set_face_handle( new_he_2, m.face_handle(he_3));
		m.set_face_handle( new_he_3, m.face_handle(he_5));

		m.set_halfedge_handle( m.face_handle(he_1), new_he_1);
		m.set_halfedge_handle( m.face_handle(he_3), new_he_2);
		m.set_halfedge_handle( m.face_handle(he_5), new_he_3);

		std::cout << "eoo" << std::endl;
//Never forget when playing with topology
		m.adjust_outgoing_halfedge( v );
		m.adjust_outgoing_halfedge( v1 );
		m.adjust_outgoing_halfedge( v2 );

		return true;
	}	

};



typedef mesh_traits<OpenMeshExtended> OpenMeshXTraits;

#include "OpenMeshX.tcc"
		  
#endif /* OPENMESHX_H_ */
