#ifdef OPENMESHX_H_
#define OPENMESHX_H_

#include <OpenMesh/Core/Utils/GenProg.hh>

template<>
class mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ox_ve_iterator : public OpenMesh::PolyMesh_ArrayKernelT<MyTraits>::VEIter
{
	public: 
                typedef OpenMesh::PolyConnectivity& mesh_ref;

                ox_ve_iterator() : 
                        OpenMesh::Iterators::VertexEdgeIterT< OpenMesh::PolyConnectivity >() {};
                
                ox_ve_iterator (mesh_ref _mesh, OpenMesh::PolyConnectivity::VertexHandle _start, bool _end) :
                        OpenMesh::Iterators::VertexEdgeIterT< OpenMesh::PolyConnectivity >(_mesh, _start, _end) {};

                ox_ve_iterator(mesh_ref _mesh, OpenMesh::PolyConnectivity::HalfedgeHandle _heh, bool _end) : 
                        OpenMesh::Iterators::VertexEdgeIterT< OpenMesh::PolyConnectivity >(_mesh, _heh, _end) {};

                ox_ve_iterator(const OpenMesh::Iterators::VertexEdgeIterT< OpenMesh::PolyConnectivity >& _rhs) :
                        OpenMesh::Iterators::VertexEdgeIterT< OpenMesh::PolyConnectivity >(_rhs) {};

                bool operator==(const ox_ve_iterator& _rhs) const 
                {
                        return 
                        ((mesh_   == _rhs.mesh_) &&
                        (start_  == _rhs.start_) &&
                        (heh_    == _rhs.heh_));
                }
 
 
                bool operator!=(const ox_ve_iterator& _rhs) const
                {
                        return !operator==(_rhs);
                }

                edge_descriptor operator*() { return this->handle(); }
                edge_descriptor operator->() { return this->handle(); }


};
                 
template<>
class mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ox_fv_iterator : public OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >
{
	public:
		typedef OpenMesh::PolyConnectivity& mesh_ref;

		ox_fv_iterator() :
			OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >() {};

		ox_fv_iterator (mesh_ref _mesh, OpenMesh::PolyConnectivity::FaceHandle _start, bool _end) :
			OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >(_mesh, _start, _end) {};

		ox_fv_iterator(mesh_ref _mesh, OpenMesh::PolyConnectivity::HalfedgeHandle _heh, bool _end) :
			OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >(_mesh, _heh, _end) {};

		ox_fv_iterator(const OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >& _rhs) :
			OpenMesh::Iterators::FaceVertexIterT< OpenMesh::PolyConnectivity >(_rhs) {};

		bool operator==(const ox_fv_iterator& _rhs) const
		{
			return
			((mesh_   == _rhs.mesh_) &&
			(start_  == _rhs.start_) &&
			(heh_    == _rhs.heh_) &&
			(lap_counter_ == _rhs.lap_counter_));
		}


		bool operator!=(const ox_fv_iterator& _rhs) const
		{
			return !operator==(_rhs);
		}

		vertex_descriptor operator*() { return this->handle(); }
		vertex_descriptor operator->() { return this->handle(); }
};
             
template <>
class mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ox_vv_iterator : public OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >
{
	public:
		typedef OpenMesh::PolyConnectivity& mesh_ref;

		ox_vv_iterator() :
			OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >() {};

		ox_vv_iterator (mesh_ref _mesh, OpenMesh::PolyConnectivity::VertexHandle _start, bool _end) :
			OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >(_mesh, _start, _end) {};

		ox_vv_iterator(mesh_ref _mesh, OpenMesh::PolyConnectivity::HalfedgeHandle _heh, bool _end) :
			OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >(_mesh, _heh, _end) {};

		ox_vv_iterator(const OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >& _rhs) :
			OpenMesh::Iterators::VertexVertexIterT< OpenMesh::PolyConnectivity >(_rhs) {};

		bool operator==(const ox_vv_iterator& _rhs) const
		{
			return
			((mesh_   == _rhs.mesh_) &&
			(start_  == _rhs.start_) &&
			(heh_    == _rhs.heh_) &&
			(lap_counter_ == _rhs.lap_counter_));
		}


		bool operator!=(const ox_vv_iterator& _rhs) const
		{
			return !operator==(_rhs);
		}

		vertex_descriptor operator*() { return this->handle(); }
		vertex_descriptor operator->() { return this->handle(); }

};
             
template <>
class mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ox_fe_iterator : public OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >
{
	public:
		typedef OpenMesh::PolyConnectivity& mesh_ref;

		ox_fe_iterator() :
			OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >() {};

		ox_fe_iterator (mesh_ref _mesh, OpenMesh::PolyConnectivity::FaceHandle _start, bool _end) :
			OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >(_mesh, _start, _end) {};

		ox_fe_iterator(mesh_ref _mesh, OpenMesh::PolyConnectivity::HalfedgeHandle _heh, bool _end) :
			OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >(_mesh, _heh, _end) {};

		ox_fe_iterator(const OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >& _rhs) :
			OpenMesh::Iterators::FaceEdgeIterT< OpenMesh::PolyConnectivity >(_rhs) {};

		bool operator==(const ox_fe_iterator& _rhs) const
		{
			return
			((mesh_   == _rhs.mesh_) &&
			(start_  == _rhs.start_) &&
			(heh_    == _rhs.heh_) &&
			(lap_counter_ == _rhs.lap_counter_));
		}


		bool operator!=(const ox_fe_iterator& _rhs) const
		{
			return !operator==(_rhs);
		}

		edge_descriptor operator*() { return this->handle(); }
		edge_descriptor operator->() { return this->handle(); }

};


template<>
std::pair
<
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::fv_iterator,
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::fv_iterator
>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_surrounding_vertices(const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_, mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_descriptor fd)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);
	return std::make_pair(m.fv_begin(fd), m.fv_end(fd));
}


template<>
std::pair
<
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::fe_iterator,
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::fe_iterator
>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_surrounding_edges(
	const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_descriptor fd)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);
	return std::make_pair(m.fe_begin(fd), m.fe_end(fd));
}


template<>
bool 
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::remove_vertex(
		mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v,
		OpenMesh::PolyMesh_ArrayKernelT<MyTraits> &m)
{
//	static_assert(TTraits::VertexAttributes != 0);
	 m.delete_vertex(v);
	 return true;
}


template<>
inline bool mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::create_face(
				  typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor a,
				  typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor b,
				  typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor c,
		  	  	  OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m)
{
	/*std::vector<mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor>  face_vhandles;


	face_vhandles.clear();
	face_vhandles.push_back(a);
	face_vhandles.push_back(b);
	face_vhandles.push_back(c);*/
	m.add_face(a, b, c);
	//m.add_face(face_vhandles);
	return true;
}

template<>
bool mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::remove_face(
				  typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_descriptor f,
		  	  	  OpenMesh::PolyMesh_ArrayKernelT<MyTraits> &m)
{
	  m.delete_face(f);
	  return true;
}

template<>
std::pair<typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_iterator,
	  	  	typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_iterator>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_all_vertices(const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_)
{
	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return std::make_pair(m.vertices_begin(), m.vertices_end());
}

template<>
std::pair<typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::edge_iterator,
	  	  	typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::edge_iterator>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_all_edges(const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_)
{
	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return std::make_pair(m.edges_begin(), m.edges_end());
}


template<>
std::pair<typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_iterator,
	  	  	typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_iterator>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_all_faces(const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_)
{
	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return std::make_pair(m.faces_begin(), m.faces_end());
}

//=========== VERTEX ADJACENCY CONCEPT ===========

template<>
bool mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::is_isolated(const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
		mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v)
{
	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return m.is_isolated(v);
}

template<>
std::pair<typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vv_iterator,
	  	  	typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vv_iterator>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_adjacent_vertices(
		const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
		  mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v)
{
	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return std::make_pair(m.vv_begin(v), m.vv_end(v));
}

template<>
std::pair<typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ve_iterator,typename mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::ve_iterator>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_adjacent_edges(
		const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
		mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v)
{

	  typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	  Mesh& m = const_cast<Mesh&>(m_);
	  return std::make_pair(m.ve_begin(v),m.ve_end(v));
}
	  
template<>
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::normal
mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_face_normal(
	const OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_descriptor f)
{
		typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
		Mesh& m = const_cast<Mesh&>(m_);
		return m.calc_face_normal(f);
}

//=========== HELPERS ==============


template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::h_edge_descriptor
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::get_previous_halfedge(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::h_edge_descriptor heh)
{
	auto t_heh = heh;

		for (t_heh = m.next_halfedge_handle(t_heh);
		m.next_halfedge_handle(t_heh) != heh;
		t_heh = m.next_halfedge_handle(t_heh) )
		{}

	return t_heh;
}


template<>
bool
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::triangulate_face(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::face_descriptor& f
)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);

	m.triangulate(f);
}


template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::calc_middle_point(
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point p1,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point p2,
	float coeff)
{
	auto diff = p2;
	diff -= p1;
	diff *= coeff;
	diff += p2;
	return diff;
}

template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::calc_translated_point(
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point p,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Point dir,
	float dist)
{
	auto diff = dir.normalize();
	diff *= dist;
	diff += p;
	return diff;
}



template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::calc_middle_vertex(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v1,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v2,
	float coeff)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);

	auto diff = m.point( v2 );
	diff -= m.point( v1 );
	diff *= coeff;
	diff += m.point( v1 );


	return m.new_vertex( diff );
}

//=========== MESH TOPOLOGY OPERATIONS ===========


template<>
bool
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::split_edge(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::edge_descriptor e)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);

	auto heh     = m.halfedge_handle(e, 0);
	auto opp_heh = m.halfedge_handle(e, 1);

	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::h_edge_descriptor new_heh, opp_new_heh, t_heh;
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor   vh;
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor   vh1( m.to_vertex_handle(heh));

	// new vertex
	vh                = advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::calc_middle_vertex
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

template<>
class advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Scalable_face
{
public:
	Scalable_face(){}
	virtual bool rescale(float coeff) = 0;
	virtual bool rescale_by_unit(float coeff) = 0;
};

template<>
class advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Scalable_face_truncate : Scalable_face
{
public:
	fv_iterator fv_begin_;
	fv_iterator fv_end_;
	Point p_;
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_;

	Scalable_face_truncate(
		OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& _m,
		std::pair<fv_iterator, fv_iterator> _fv_pair,
		Point _p) : m_( _m )
	{
		p_ = _p;
		fv_begin_ = _fv_pair.first;
		fv_end_ = _fv_pair.second;
	}

	bool rescale(float coeff)
	{
		for (auto fv_i = fv_begin_; fv_i != fv_end_; ++fv_i)
		{
			auto p_new = calc_middle_point(p_, m_.point(*fv_i), coeff);
			m_.set_point(*fv_i, p_new);
		}
	}

	bool rescale_by_unit(float coeff)
	{
		for (auto fv_i = fv_begin_; fv_i != fv_end_; ++fv_i)
		{
			auto dir= m_.point(*fv_i);
			dir -= p_;
			auto p_new = calc_translated_point(m_.point(*fv_i), dir, coeff);
			m_.set_point(*fv_i, p_new);
		}

	}
};


template<>
class advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Scalable_face_bevel : Scalable_face
{
public:
	Point p1_;
	Point p2_;
	std::vector<vertex_descriptor> vector1_;
	std::vector<vertex_descriptor> vector2_;
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_;

	Scalable_face_bevel(
		OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& _m,
		Point _p1,
		Point _p2,
		std::vector<vertex_descriptor> _vector1,
		std::vector<vertex_descriptor> _vector2
		) :
			m_( _m ),
			vector1_( _vector1 ),
			vector2_( _vector2 )
	{
		p1_ = _p1;
		p2_ = _p2;
	}

	bool rescale(float coeff)
	{
		for (auto vertex : vector1_)
		{
			auto p_new = calc_middle_point(p1_, m_.point(vertex), coeff);
			m_.set_point(vertex, p_new);
		}
		for (auto vertex : vector2_)
		{
			auto p_new = calc_middle_point(p2_, m_.point(vertex), coeff);
			m_.set_point(vertex, p_new);
		}
	}

	bool rescale_by_unit(float coeff)
	{
		for (auto vertex : vector1_)
		{
			auto dir = m_.point( vertex );
			dir -= p1_;
			auto p_new = calc_translated_point(m_.point(vertex), dir, coeff);
			m_.set_point(vertex, p_new);

		}
		for (auto vertex : vector2_)
		{
			auto dir = m_.point( vertex );
			dir -= p2_;
			auto p_new = calc_translated_point(m_.point( vertex ), dir, coeff);
			m_.set_point(vertex, p_new);

		}

	}
};

template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Scalable_face_truncate
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::truncate(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::vertex_descriptor v,
	float coeff)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
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

	auto res_p = m.point(v);

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

	return Scalable_face_truncate(m_, get_surrounding_vertices(m, f_a), res_p);
}	


template<>
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::Scalable_face_bevel
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::bevel(
		OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
		advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::edge_descriptor e,
		float coeff)
{
	
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
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

	auto pv_a = m.point( v_a );
	auto pv_b = m.point( v_b );

	auto prev_he_1 = he_a;
	auto he_1 = m.next_halfedge_handle(prev_he_1);
	auto opp_he_1 = m.opposite_halfedge_handle(he_1);
	auto v_new = v_a;

	std::vector<vertex_descriptor> vd_list_a;
	std::vector<vertex_descriptor> vd_list_b;

	r_l.push_back(std::make_tuple(v_new, prev_he_1, he_1, opp_he_1));
	vd_list_a.push_back(v_new);

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
		vd_list_a.push_back(v_new);
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

	vd_list_b.push_back(v_new);

	do
	{	
		prev_he_1 = opp_he_1;
		he_1 = m.next_halfedge_handle(opp_he_1);
		opp_he_1 = m.opposite_halfedge_handle(he_1);

		if (m.next_halfedge_handle(opp_he_1) != he_a)
		{
			v_new = calc_middle_vertex(
				m,
				v_b,
				m.to_vertex_handle(he_1),
				coeff);
			vd_list_b.push_back(v_new);
		}
		else
		{
			v_new = v_b;
			vd_list_b.push_back(v_new);
		}

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

	return Scalable_face_bevel(m, pv_a, pv_b,
		vd_list_a, vd_list_b);
}


template<>
bool
advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::fill_ring(
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits>& m_,
	advanced_mesh_traits<OpenMesh::PolyMesh_ArrayKernelT<MyTraits>>::edge_descriptor e)
{
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
	Mesh& m = const_cast<Mesh&>(m_);

	std::vector<h_edge_descriptor> inside_face_he;

	auto he_a = m.halfedge_handle(e, 0);
	auto he_b = m.halfedge_handle(e, 1);
	if (m.is_boundary(he_b))
		std::swap(he_a, he_b);

	auto f_a = m.new_face();

	auto he_n = get_previous_halfedge(m, he_b);
	auto he_o = m.opposite_halfedge_handle(he_n);
	auto he_m = get_previous_halfedge(m, he_o);
	auto he_next = m.opposite_halfedge_handle(he_m);

		while (!m.is_boundary(he_next))
		{
			he_m = get_previous_halfedge(m, he_next);
			he_next = m.opposite_halfedge_handle(he_m);
		}

	inside_face_he.push_back(he_next);
	m.set_face_handle(he_next, f_a);
	m.set_halfedge_handle(f_a, he_next);

	do
	{
		he_n = get_previous_halfedge(m, he_m);
		he_o = m.opposite_halfedge_handle(he_n);
		he_m = get_previous_halfedge(m, he_o);
		he_next = m.opposite_halfedge_handle(he_m);

		while (!m.is_boundary(he_next))
		{
			he_m = get_previous_halfedge(m, he_next);
			he_next = m.opposite_halfedge_handle(he_m);
		}

		inside_face_he.push_back(he_next);
		m.set_face_handle(he_next, f_a);
	}
	while (he_next != he_a);

	for (int i=0; i<inside_face_he.size(); i++)
	{
		auto hx_a = inside_face_he[ i ];
		auto hx_b = inside_face_he[ (i+1) % inside_face_he.size()];
		m.set_next_halfedge_handle(hx_a, hx_b);
	}

	return true;
}



#endif // OPENMESHX_H_
