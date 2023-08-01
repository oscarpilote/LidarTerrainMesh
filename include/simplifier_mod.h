#include <stddef.h>

/* Experimental: Mesh simplifier
 * Reduces the number of triangles in the mesh, attempting to preserve mesh
 * appearance as much as possible The algorithm tries to preserve mesh topology
 * and can stop short of the target goal based on topology constraints or target
 * error. If not all attributes from the input mesh are required, it's
 * recommended to reindex the mesh using meshopt_generateShadowIndexBuffer prior
 * to simplification. Returns the number of indices after simplification, with
 * destination containing new index data The resulting index buffer references
 * vertices from the original vertex buffer. If the original vertex data isn't
 * required, creating a compact vertex buffer using meshopt_optimizeVertexFetch
 * is recommended.
 *
 * destination must contain enough space for the target index buffer, worst case
 * is index_count elements (*not* target_index_count)! vertex_positions should
 * have float3 position in the first 12 bytes of each vertex - similar to
 * glVertexPointer target_error represents the error relative to mesh extents
 * that can be tolerated, e.g. 0.01 = 1% deformation result_error can be NULL;
 * when it's not NULL, it will contain the resulting (relative) error after
 * simplification The optional simplificationn_remap allows to keep track of the
 * collapses. If not NULL it must contain sufficient space for vertex_count
 * unsigned int.
 */
size_t meshopt_simplify_mod(unsigned int *dest_idx, unsigned int *simp_remap,
			    const unsigned int *indices, size_t index_count,
			    const float *vertex_positions, size_t vertex_count,
			    size_t vertex_positions_stride,
			    size_t target_index_count, float target_error,
			    float *result_error);

