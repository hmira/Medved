#ifndef __RANGE_EDGES_HPP__
#define __RANGE_EDGES_HPP__

#include <type_traits>
#include <hmira/range/range.hpp>

namespace hmira
{

namespace range
{

	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	TFunction
	for_each_edge(TMesh& m, TFunction a_unary_function)
	{
		auto all_edges = TMesh_Traits::get_all_edges(m);
		for (auto f = all_edges.first; f != all_edges.second; ++f)
		{
			a_unary_function(f);
		}
		return a_unary_function;
	}
	
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	typename TMesh_Traits::edge_iterator
	find_edge_if(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nfind_edge_if() requires functor that returns bool\n******\n");
		
		auto all_edges = TMesh_Traits::get_all_edges(m);
		auto edge = all_edges.first;
		while ( edge != all_edges.second)
		{
			if (a_unary_predicate(*edge)) return edge;
			++edge;
		}
		return all_edges.second;
	}
	
		
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	bool
	all_of_edges(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nall_of_edges() requires functor that returns bool\n******\n");

		auto all_edges = TMesh_Traits::get_all_edges(m);
		auto edge = all_edges.first;
		
		while ( edge != all_edges.second)
		{
			if (!a_unary_predicate(*edge)) return false;
			++edge;
		}
		return true;
	}
	
	
} // namespace range
	
} // namespace hmira

#endif // __RANGE_EDGES_HPP__