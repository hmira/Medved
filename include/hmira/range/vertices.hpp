#ifndef __RANGE_VERTICES_HPP__
#define __RANGE_VERTICES_HPP__

#include <type_traits>
#include <hmira/range/range.hpp>

namespace hmira
{

namespace range
{

	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	TFunction
	for_each_vertex(TMesh& m, TFunction a_unary_function)
	{
		auto all_vertices = TMesh_Traits::get_all_vertices(m);
		for (auto f = all_vertices.first; f != all_vertices.second; ++f)
		{
			a_unary_function(f);
		}
		return a_unary_function;
	}
	
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	typename TMesh_Traits::vertex_iterator
	find_vertex_if(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nfind_vertex_if() requires functor that returns bool\n******\n");
		
		auto all_vertices = TMesh_Traits::get_all_vertices(m);
		auto vertex = all_vertices.first;
		while ( vertex != all_vertices.second)
		{
			if (a_unary_predicate(*vertex)) return vertex;
			++vertex;
		}
		return all_vertices.second;
	}
	
		
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	bool
	all_of_vertices(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nall_of_vertices() requires functor that returns bool\n******\n");

		auto all_vertices = TMesh_Traits::get_all_vertices(m);
		auto vertex = all_vertices.first;
		
		while ( vertex != all_vertices.second)
		{
			if (!a_unary_predicate(*vertex)) return false;
			++vertex;
		}
		return true;
	}
	
	
} // namespace range
	
} // namespace hmira

#endif // __RANGE_VERTICES_HPP__