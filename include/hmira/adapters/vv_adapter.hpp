#ifndef __VV_ADAPTER_HPP__
#define __VV_ADAPTER_HPP__

#include <deque>

#include <hmira/range/faces.hpp>

#include <hmira/preprocessor/debug_mode.hpp>

namespace hmira
{

namespace adapters
{
	template <typename TMesh, typename TMesh_Traits>
	class vv_adapter
	{
	public:
		typedef typename std::deque<typename TMesh_Traits::vertex_descriptor>::iterator vv_iterator;
		static inline std::pair<
			vv_iterator,
			vv_iterator>
		get_adjacent_vertices(
				const TMesh& m_,
				typename TMesh_Traits::vertex_descriptor v)
				{
					return vv_iterator_adapter(m_, v);
				}

		static std::pair<
			typename std::deque<
					typename TMesh_Traits::vertex_descriptor
				>::iterator,
				typename std::deque<
					typename TMesh_Traits::vertex_descriptor
				>::iterator
			>
		vv_iterator_adapter(const TMesh& m, typename TMesh_Traits::vertex_descriptor tested_vertex)
		{
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : ---------------------------------" )
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : acquiring vv_iterator started" )
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : ---------------------------------" )
			
			typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
			typedef typename TMesh_Traits::face_descriptor face_descriptor;
			std::deque<vertex_descriptor> vertex_vector;
			std::deque<vertex_descriptor> vertices_of_tested_face;
		
			auto all_faces = TMesh_Traits::get_all_faces(m);
			for (auto f = all_faces.first; f != all_faces.second; ++f)
			{
				auto face = *f;
				int found = 0;
				auto sur_v_face = TMesh_Traits::get_surrounding_vertices(m, face);

				for (auto vi = sur_v_face.first; vi != sur_v_face.second; ++vi)
				{
					if (*vi == tested_vertex)
					{
						auto found_v = vi;
						++found_v;
						if (found_v !=	sur_v_face.second)
						{
							vertex_vector.push_back(*found_v);
						}
						else
						{
							vertex_vector.push_back(*sur_v_face.first);
						}
					}
				}

			}
			
			auto first = vertex_vector.begin();
			auto second = vertex_vector.end();
			
			auto result = std::make_pair(first, second);
			
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : ---------------------------------" )
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : the pair of vv_iterators acquired" )
			H_DEBUG_STDERR( "[VERTEX-VERTEX ITERATOR ADAPTER] : ---------------------------------" )
			
			return result;
		}
		
		
	};

}//adapters

}//hmira

#endif // __VV_ADAPTER_HPP__