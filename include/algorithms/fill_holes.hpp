
template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int fill_holes(TMesh& m)
{
	typedef typename TMesh_Traits::edge_descriptor edge_descriptor;
        typedef typename TMesh_Traits::edge_iterator edge_iterator;

        auto all_edges = TMesh_Traits::get_all_edges(m);

	for (auto i = all_edges.first; i != all_edges.second; ++i)
	{
		edge_descriptor ed = *i;
		if (TMesh_Traits::is_boundary(m, ed)) TMesh_Traits::fill_ring(m, ed);
	}

	return 0;
}