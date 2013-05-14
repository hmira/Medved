#ifndef __TRIANGULATE_HPP__
#define __TRIANGULATE_HPP__

template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int triangulate(TMesh& m)
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
        typedef typename TMesh_Traits::face_iterator face_iterator;
	typedef typename TMesh_Traits::vertex_descriptor vertex_descriptor;
	typedef typename TMesh_Traits::fv_iterator fv_iterator;

        auto all_faces = TMesh_Traits::get_all_faces(m);

	for (auto i = all_faces.first; i != all_faces.second; ++i)
	{
		face_descriptor fd = *i;
		TMesh_Traits::triangulate_face(m, fd);
	}

	return 0;
}

#endif //__TRIANGULATE_HPP__