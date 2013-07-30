#ifndef __FF_ADAPTER_HPP__
#define __FF_ADAPTER_HPP__

#include <deque>

#include <hmira/range/faces.hpp>

#include <hmira/preprocessor/debug_mode.hpp>

namespace hmira
{

namespace adapters
{
	template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
	class ff_adapter
	{
	public:
		
		typedef typename std::deque<typename TMesh_Traits::face_descriptor>::iterator ffiterator;
		static
		std::pair<
			typename std::deque<
					typename TMesh_Traits::face_descriptor
				>::iterator,
				typename std::deque<
					typename TMesh_Traits::face_descriptor
				>::iterator
			>
		ff_iterator_adapter(TMesh& m, typename TMesh_Traits::face_descriptor tested_face)
		{
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : ---------------------------------" )
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : acquiring ff_iterator started" )
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : ---------------------------------" )
			
			typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
			typedef typename TMesh_Traits::face_descriptor face_descriptor;
			std::deque<face_descriptor> face_vector;
			std::deque<vertex_descriptor> vertices_of_tested_face;
			
			auto sur_v_tested = TMesh_Traits::get_surrounding_vertices(m, tested_face);
			for (auto vi = sur_v_tested.first; vi != sur_v_tested.second; ++vi)
			{
				vertices_of_tested_face.push_back(*vi);
			}
			
			range::for_each_face(m,
				[&]
				(typename TMesh_Traits::face_descriptor face)
				{
					int found = 0;
					auto sur_v_face = TMesh_Traits::get_surrounding_vertices(m, face);
					for (auto vt : vertices_of_tested_face)
					{
						for (auto vi = sur_v_face.first; vi != sur_v_face.second; ++vi)
						{
							if (*vi == vt)
							{
								found++;
							}
						}
					}
					if (found == 2) //if 2 vertices are in common, the faces are in neighbourhood
					{
						face_vector.push_back(face);
					}
				});
			
			auto first = face_vector.begin();
			auto second = face_vector.end();
			
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : ---------------------------------" )
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : the pair of ff_iterators acquired" )
			H_DEBUG_STDERR( "[FACE-FACE ITERATOR ADAPTER] : ---------------------------------" )
			
			auto result = std::make_pair(first, second);
			
			return result;
		}
	};
}//adapters

}//hmira

#endif // __FF_ADAPTER_HPP__