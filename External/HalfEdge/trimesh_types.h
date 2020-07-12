#ifndef __trimesh_types_h__
#define __trimesh_types_h__

namespace trimesh
{
    typedef long index_t;
    
    struct edge_t
    {
        index_t v[2];
        
        index_t& start() { return v[0]; }
        const index_t& start() const { return v[0]; }
        
        index_t& end() { return v[1]; }
        const index_t& end() const { return v[1]; }
        
        edge_t()
        {
            v[0] = v[1] = -1;
        }
    };
    
    struct triangle_t
    {
		union
		{
			struct
			{
				index_t i, j, k;
			};
			index_t v[3];
		};
        
		index_t& operator[](int index) { return v[index]; }
		const index_t& operator[](int index) const { return v[index]; }

        triangle_t()
        {
            v[0] = v[1] = v[2] = -1;
        }
    };
}

#endif /* __trimesh_types_h__ */
