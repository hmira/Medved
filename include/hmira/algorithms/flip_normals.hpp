#ifndef __FLIP_NORMALS_HPP__
#define __FLIP_NORMALS_HPP__

template <typename TMesh, typename TMesh_Traits = mesh_traits<TMesh>>
int flip_normals(TMesh& m)
{
	typedef typename TMesh_Traits::face_descriptor face_descriptor;
	typedef typename TMesh_Traits::face_iterator face_iterator;

	auto all_faces = TMesh_Traits::get_all_faces(m);

	for (auto i = all_faces.first; i != all_faces.second; ++i)
	{
		face_descriptor fd = *i;
		TMesh_Traits::flip_face_normal(m, fd);
	}

	return 0;
}

#endif //__FLIP_NORMALS_HPP__