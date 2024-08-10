#include "base_mesh.h"

namespace SPH
{
//=============================================================================================//
Arrayi Mesh::transfer1DtoMeshIndex(const Arrayi &mesh_size, size_t i) const
{
    size_t row_size = mesh_size[1];
    size_t column = i / row_size;
    return Arrayi(column, i - column * row_size);
}
//=============================================================================================//
size_t Mesh::transferMeshIndexToMortonOrder(const Arrayi &mesh_index) const
{
    return MortonCode(mesh_index[0]) | (MortonCode(mesh_index[1]) << 1);
}
} // namespace SPH
//=============================================================================================//
