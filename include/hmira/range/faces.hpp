#ifndef __RANGE_FACES_HPP__
#define __RANGE_FACES_HPP__

#include <type_traits>
#include <hmira/range/range.hpp>

namespace hmira
{

namespace range
{

	/**
	 * @param m	mesh
	 * @param a_unary_function function called after every iteration
	 */
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	TFunction
	for_each_face(TMesh& m, TFunction a_unary_function)
	{
		auto all_faces = TMesh_Traits::get_all_faces(m);
		for (auto f = all_faces.first; f != all_faces.second; ++f)
		{
			a_unary_function(f);
		}
		return a_unary_function;
	}
	
	/**
	 * @param m	mesh
	 * @param a_unary_function boolean function
	 * @return found iterator, if nothing found, returns last iterator
	 */
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	typename TMesh_Traits::face_iterator
	find_face_if(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nfind_face_if() requires functor that returns bool\n******\n");
		
		auto all_faces = TMesh_Traits::get_all_faces(m);
		auto face = all_faces.first;
		while ( face != all_faces.second)
		{
			if (a_unary_predicate(*face)) return face;
			++face;
		}
		return all_faces.second;
	}
	
		
	template <typename TFunction, typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	inline
	bool
	all_of_faces(TMesh& m, TFunction a_unary_predicate)
	{
		typedef function_traits<decltype(a_unary_predicate)> traits;
		typedef typename traits::result_type lambda_type;
		static_assert(
			std::is_same<bool, lambda_type>::value,
			"\n******\nall_of_faces() requires functor that returns bool\n******\n");

		auto all_faces = TMesh_Traits::get_all_faces(m);
		auto face = all_faces.first;
		
		while ( face != all_faces.second)
		{
			if (!a_unary_predicate(*face)) return false;
			++face;
		}
		return true;
	}
	
	
} // namespace range
	
} // namespace hmira

#endif // __RANGE_FACES_HPP__