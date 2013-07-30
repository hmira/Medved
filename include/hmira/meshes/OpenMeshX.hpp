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
#include <algorithm>
#include <memory>
#include <functional>

#include <boost/mpl/if.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/int.hpp>

#include <hmira/adapters/vv_adapter.hpp>

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
template<typename TTraits>
class mesh_traits< OpenMesh::PolyMesh_ArrayKernelT<TTraits> /*OpenMeshExtended*/>

{
public:
	typedef OpenMesh::PolyMesh_ArrayKernelT<TTraits> OpenMeshExtended;

	typedef OpenMesh::DefaultTraits original_traits;

	typedef typename OpenMeshExtended::Point Point;
	typedef typename OpenMeshExtended::Normal normal;


	typedef typename OpenMeshExtended::HalfedgeHandle h_edge_descriptor;//consider moving to advanced traits
	typedef typename OpenMeshExtended::HalfedgeIter h_edge_iterator;
	//typedef typename std::vector<OpenMeshExtended::Halfedge>::size_type h_edges_size_type;

	typedef typename OpenMeshExtended::VertexHandle vertex_descriptor;
	typedef typename OpenMeshExtended::VertexIter vertex_iterator;
	//typedef typename std::vector<OpenMeshExtended::Vertex>::size_type vertices_size_type;

	typedef typename OpenMeshExtended::EdgeHandle edge_descriptor;
	typedef typename OpenMeshExtended::EdgeIter edge_iterator;
	//typedef typename std::vector<OpenMeshExtended::Edge>::size_type edges_size_type;

	typedef typename OpenMeshExtended::FaceHandle face_descriptor;
	typedef typename OpenMeshExtended::FaceIter face_iterator;
	//typedef typename std::vector<OpenMeshExtended::Face>::size_type faces_size_type;

	class ox_fv_iterator;
	class ox_vv_iterator;
	class ox_ve_iterator;
	class ox_fe_iterator;

	typedef ox_fv_iterator fv_iterator;
// 	typedef ox_vv_iterator vv_iterator;
// 	typedef typename hmira::adapters::vv_adapter<OpenMeshExtended, mesh_traits<OpenMeshExtended> >::vv_iterator vv_iterator;

	typedef ox_ve_iterator ve_iterator;
	typedef ox_fe_iterator fe_iterator;

	static std::pair<fv_iterator, fv_iterator>
	get_surrounding_vertices(const OpenMeshExtended& m_, face_descriptor fd);

	static std::pair<fe_iterator, fe_iterator>
	get_surrounding_edges(const OpenMeshExtended& m_, face_descriptor fd);

//=================CONCEPTS======================

static vertex_descriptor create_vertex(
				Point p,
				OpenMeshExtended &m) {return m.add_vertex(p);}
static bool remove_vertex(
				vertex_descriptor v,
		  	  	OpenMeshExtended &m);

inline static bool create_face(
				  vertex_descriptor a,
				  vertex_descriptor b,
				  vertex_descriptor c,
		  	  	  OpenMeshExtended& m);

static bool remove_face(
				  typename mesh_traits<OpenMeshExtended>::face_descriptor f,
		  	  	  OpenMeshExtended& m);

static inline
std::pair<vertex_descriptor, vertex_descriptor>
	get_edge_vertices(const OpenMeshExtended& m_, edge_descriptor e)
	{
		auto he0 = m_.halfedge_handle(e, 0);
		auto he1 = m_.halfedge_handle(e, 1);
		auto v0 = m_.to_vertex_handle(he0);
		auto v1 = m_.to_vertex_handle(he1);
		return std::make_pair<vertex_descriptor, vertex_descriptor>(
			std::move(v0),
			std::move(v1));
	}

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

// static std::pair<vv_iterator,
// 	  	  vv_iterator>
// get_adjacent_vertices(
// 		const OpenMeshExtended& m_,
// 		  vertex_descriptor v){}

static std::pair<vv_iterator,
		vv_iterator>
get_adjacent_vertices(
	const OpenMeshExtended& m_,
		vertex_descriptor v)
{
	return hmira::adapters::vv_adapter<OpenMeshExtended, mesh_traits<OpenMeshExtended> >::vv_iterator_adapter(m_, v);
}

static std::pair<ve_iterator, ve_iterator>
get_adjacent_edges(
		const OpenMeshExtended& m_,
		vertex_descriptor v);

//============ NORMAL ===================
 
static typename mesh_traits<OpenMeshExtended>::normal 
get_face_normal(
	const OpenMeshExtended& m_,
	face_descriptor f);

//============ PROPERTIES ==================
template <typename TProperty>
static
void
set_property(
	const OpenMeshExtended& m_,
	vertex_descriptor v,
	TProperty p){}

//=============== COORDINATES =============
static
inline Point
get_coordinates(
	const OpenMeshExtended& m_,
	vertex_descriptor v)	{ return m_.point(v); }

};

template<typename TTraits>
class advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<TTraits>> : public mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<TTraits>>
{
public:

	typedef OpenMesh::PolyMesh_ArrayKernelT<TTraits> OpenMeshExtended;

	typedef OpenMesh::DefaultTraits original_traits;
	typedef mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<TTraits>> parent_traits;

	typedef typename parent_traits::Point Point;
	typedef typename parent_traits::normal normal;


	typedef typename parent_traits::h_edge_descriptor h_edge_descriptor;
	typedef typename parent_traits::h_edge_iterator h_edge_iterator;

	typedef typename parent_traits::vertex_descriptor vertex_descriptor;
	typedef typename parent_traits::vertex_iterator vertex_iterator;

	typedef typename parent_traits::edge_descriptor edge_descriptor;
	typedef typename parent_traits::edge_iterator edge_iterator;

	typedef typename parent_traits::face_descriptor face_descriptor;
	typedef typename parent_traits::face_iterator face_iterator;

	typedef typename parent_traits::fv_iterator fv_iterator;
	typedef typename parent_traits::vv_iterator vv_iterator;
// 	typedef typename hmira::adaptors::vv_iterator_adaptor_traits<OpenMeshExtended, mesh_traits<OpenMeshExtended> >::vv_iterator vv_iterator;
	typedef typename parent_traits::ve_iterator ve_iterator;


	
	struct Face_has_normal_tag {};
	struct Face_has_no_normal_tag {};
	typedef typename boost::mpl::if_<
		boost::mpl::bitand_<
			boost::mpl::int_<TTraits::FaceAttributes>,
			boost::mpl::int_<OpenMesh::Attributes::Normal>
		>, 
		Face_has_normal_tag, 
		Face_has_no_normal_tag
	>::type Face_normal_tag;

	static bool 
	flip_face_normal_t(
		OpenMeshExtended& m_,
		face_descriptor f,
		Face_has_normal_tag
  			)
        {

		typedef OpenMeshExtended Mesh;
		Mesh& m = const_cast<Mesh&>(m_);

		auto n = m.normal(f);
		m.set_normal(f, -n);
		return true;
	}

	static bool 
	flip_face_normal_t(
		OpenMeshExtended& m_,
		face_descriptor f,
		Face_has_no_normal_tag
  			)
        {
		std::cerr << "NO NORMAL ATTRIBUTES SET" << std::endl;
		return false;
	}

	static bool
	flip_face_normal(
		OpenMeshExtended& m_,
		face_descriptor f)
	{
		Face_normal_tag tag;
		return flip_face_normal_t(m_, f, tag);
	}

	static bool
	triangulate_face(
		OpenMeshExtended& m_,
		face_descriptor& f);


	class Scalable_face;
	class Scalable_face_truncate;
	class Scalable_face_bevel;

	static Scalable_face_truncate
	truncate(
		OpenMeshExtended& m_,
		vertex_descriptor v,
		float coeff = 0.5f);


	static Scalable_face_bevel
	bevel(
		OpenMeshExtended& m_,
		edge_descriptor e,
		  float coeff = 0.5f);

	static bool
	fill_ring(
		OpenMeshExtended& m_,
		edge_descriptor e);

//=========== HELPERS ==============
	
	static Point
	calc_middle_point(
		Point p1,
		Point p2,
		float coeff = 0.5f);

	static Point
	calc_translated_point(
		Point p_from,
		Point p_to,
		float dist = 0.5f);


	static vertex_descriptor
	calc_middle_vertex(
		OpenMeshExtended& m_,
		vertex_descriptor v1,
		vertex_descriptor v2,
		float coeff = 0.5f);

	static bool
	split_edge(
		OpenMeshExtended& m_,
		edge_descriptor e);
	
	static bool
	is_boundary(
		OpenMeshExtended& m_,
		edge_descriptor e)	{return m_.is_boundary(e);}

	static
	h_edge_descriptor
	get_previous_halfedge(
		OpenMeshExtended& m,
		h_edge_descriptor heh);

};




typedef mesh_traits<OpenMeshExtended> OpenMeshXTraits;

#include "OpenMeshX.tcc"
		  
#endif /* OPENMESHX_H_ */
