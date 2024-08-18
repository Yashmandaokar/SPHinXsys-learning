#include "level_set.h"

#include "base_body.h"
#include "base_kernel.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "mesh_iterators.hpp"
#include "tbb/parallel_sort.h"

namespace SPH
{
//=================================================================================================//
void LevelSet::initializeDataForSingularPackage(const size_t package_index, Real far_field_level_set)
{
    auto &phi = phi_.DataField()[package_index];
    auto &near_interface_id = near_interface_id_.DataField()[package_index];
    auto &phi_gradient = phi_gradient_.DataField()[package_index];
    auto &kernel_weight = kernel_weight_.DataField()[package_index];
    auto &kernel_gradient = kernel_gradient_.DataField()[package_index];

    for_each_cell_data(
        [&](int i, int j, int k)
        {
            phi[i][j][k] = far_field_level_set;
            near_interface_id[i][j][k] = far_field_level_set < 0.0 ? -2 : 2;
            phi_gradient[i][j][k] = Vecd::Ones();
            kernel_weight[i][j][k] = far_field_level_set < 0.0 ? 0 : 1.0;
            kernel_gradient[i][j][k] = Vec3d::Zero();
        });
}
//=================================================================================================//
bool LevelSet::isWithinCorePackage(Vecd position)
{
    Arrayi cell_index = CellIndexFromPosition(position);
    return isCoreDataPackage(cell_index);
}
//=============================================================================================//
void LevelSet::diffuseLevelSetSign()
{
    diffuse_level_set_sign.exec();
    // package_parallel_for(
    //     [&](size_t package_index)
    //     {
    //         auto phi_data = phi_.DataField();
    //         auto near_interface_id_data = near_interface_id_.DataField();
    //         auto &neighborhood = cell_neighborhood_[package_index];

    //         for_each_cell_data(
    //             [&](int i, int j, int k)
    //             {
    //                 // near interface cells are not considered
    //                 if (abs(near_interface_id_data[package_index][i][j][k]) > 1)
    //                 {
    //                     mesh_find_if3d<-1, 2>(
    //                         [&](int l, int m, int n) -> bool
    //                         {
    //                             NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + l, j + m, k + n), neighborhood);
    //                             int near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
    //                             bool is_found = abs(near_interface_id) == 1;
    //                             if (is_found)
    //                             {
    //                                 Real phi_0 = phi_data[package_index][i][j][k];
    //                                 near_interface_id_data[package_index][i][j][k] = near_interface_id;
    //                                 phi_data[package_index][i][j][k] = near_interface_id == 1 ? fabs(phi_0) : -fabs(phi_0);
    //                             }
    //                             return is_found;
    //                         });
    //                 }
    //             });
    //     });
}
//=============================================================================================//
void LevelSet::reinitializeLevelSet()
{
    reinitialize_level_set.exec();
    // package_parallel_for(
    //     [&](size_t package_index)
    //     {
    //         auto phi_data = phi_.DataField();
    //         auto &phi_addrs = phi_data[package_index];
    //         auto &near_interface_id_addrs = near_interface_id_.DataField()[package_index];
    //         auto &neighborhood = cell_neighborhood_[package_index];

    //         for_each_cell_data(
    //             [&](int i, int j, int k)
    //             {
    //                 // only reinitialize non cut cells
    //                 if (near_interface_id_addrs[i][j][k] != 0)
    //                 {
    //                     Real phi_0 = phi_addrs[i][j][k];
    //                     Real sign = phi_0 / sqrt(phi_0 * phi_0 + data_spacing_ * data_spacing_);
    //                     NeighbourIndex x1 = NeighbourIndexShift(Arrayi(i + 1, j, k), neighborhood);
    //                     NeighbourIndex x2 = NeighbourIndexShift(Arrayi(i - 1, j, k), neighborhood);
    //                     NeighbourIndex y1 = NeighbourIndexShift(Arrayi(i, j + 1, k), neighborhood);
    //                     NeighbourIndex y2 = NeighbourIndexShift(Arrayi(i, j - 1, k), neighborhood);
    //                     NeighbourIndex z1 = NeighbourIndexShift(Arrayi(i, j, k + 1), neighborhood);
    //                     NeighbourIndex z2 = NeighbourIndexShift(Arrayi(i, j, k - 1), neighborhood);
    //                     Real dv_x = upwindDifference(sign, phi_data[x1.first][x1.second[0]][x1.second[1]][x1.second[2]] - phi_0, phi_0 - phi_data[x2.first][x2.second[0]][x2.second[1]][x2.second[2]]);
    //                     Real dv_y = upwindDifference(sign, phi_data[y1.first][y1.second[0]][y1.second[1]][y1.second[2]] - phi_0, phi_0 - phi_data[y2.first][y2.second[0]][y2.second[1]][y2.second[2]]);
    //                     Real dv_z = upwindDifference(sign, phi_data[z1.first][z1.second[0]][z1.second[1]][z1.second[2]] - phi_0, phi_0 - phi_data[z2.first][z2.second[0]][z2.second[1]][z2.second[2]]);
    //                     phi_addrs[i][j][k] -= 0.3 * sign * (Vec3d(dv_x, dv_y, dv_z).norm() - data_spacing_);
    //                 }
    //             });
    //     });
}
//=================================================================================================//
void LevelSet::markNearInterface(Real small_shift_factor)
{
    // Real small_shift = small_shift_factor * data_spacing_;
    mark_near_interface.setSmallShiftFactor(small_shift_factor);
    mark_near_interface.exec();

    // package_parallel_for(
    //     [&](size_t package_index)
    //     {
    //         auto &phi_addrs = phi_.DataField()[package_index];
    //         auto &near_interface_id_addrs = near_interface_id_.DataField()[package_index];
    //         auto neighborhood = cell_neighborhood_[package_index];

    //         // corner averages, note that the first row and first column are not used
    //         PackageTemporaryData<Real> corner_averages;
    //         mesh_for_each3d<0, pkg_size + 1>(
    //             [&](int i, int j, int k)
    //             {
    //                 corner_averages[i][j][k] = CornerAverage(phi_, Arrayi(i, j, k), Arrayi(-1, -1, -1), neighborhood);
    //             });

    //         for_each_cell_data(
    //             [&](int i, int j, int k)
    //             {
    //                 // first assume far cells
    //                 Real phi_0 = phi_addrs[i][j][k];
    //                 int near_interface_id = phi_0 > 0.0 ? 2 : -2;
    //                 if (fabs(phi_0) < small_shift)
    //                 {
    //                     near_interface_id = 0;
    //                     Real phi_average_0 = corner_averages[i][j][k];
    //                     // find inner and outer cut cells
    //                     mesh_for_each3d<0, 2>(
    //                         [&](int l, int m, int n)
    //                         {
    //                             Real phi_average = corner_averages[i + l][j + m][k + n];
    //                             if ((phi_average_0 - small_shift) * (phi_average - small_shift) < 0.0)
    //                                 near_interface_id = 1;
    //                             if ((phi_average_0 + small_shift) * (phi_average + small_shift) < 0.0)
    //                                 near_interface_id = -1;
    //                         });
    //                     // find zero cut cells
    //                     mesh_for_each3d<0, 2>(
    //                         [&](int l, int m, int n)
    //                         {
    //                             Real phi_average = corner_averages[i + l][j + m][k + n];
    //                             if (phi_average_0 * phi_average < 0.0)
    //                                 near_interface_id = 0;
    //                         });
    //                 }
    //                 // assign this is to package
    //                 near_interface_id_addrs[i][j][k] = near_interface_id;
    //             });
    //     });
}
//=================================================================================================//
void LevelSet::redistanceInterfaceForAPackage(const size_t package_index)
{
    auto phi_data = phi_.DataField();
    auto near_interface_id_data = near_interface_id_.DataField();
    auto &neighborhood = cell_neighborhood_[package_index];

    for_each_cell_data(
        [&](int i, int j, int k)
        {
            int near_interface_id = near_interface_id_data[package_index][i][j][k];
            if (near_interface_id == 0)
            {
                bool positive_band = false;
                bool negative_band = false;
                mesh_for_each3d<-1, 2>(
                    [&](int r, int s, int t)
                    {
                        NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + r, j + s, k + t), neighborhood);
                        int neighbor_near_interface_id = near_interface_id_data[neighbour_index.first][neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                        if (neighbor_near_interface_id >= 1)
                            positive_band = true;
                        if (neighbor_near_interface_id <= -1)
                            negative_band = true;
                    });
                if (positive_band == false)
                {
                    Real min_distance_p = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + x, j + y, k + z), neighborhood);
                            auto &neighbor_phi = phi_.DataField()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.DataField()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.DataField()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] >= 1)
                            {
                                Real phi_p_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_p = SMIN(min_distance_p, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ + phi_p_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j][k] = -min_distance_p;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j][k] = -1;
                }
                if (negative_band == false)
                {
                    Real min_distance_n = 5.0 * data_spacing_;
                    mesh_for_each3d<-4, 5>(
                        [&](int x, int y, int z)
                        {
                            NeighbourIndex neighbour_index = NeighbourIndexShift(Arrayi(i + x, j + y, k + z), neighborhood);
                            auto &neighbor_phi = phi_.DataField()[neighbour_index.first];
                            auto &neighbor_phi_gradient = phi_gradient_.DataField()[neighbour_index.first];
                            auto &neighbor_near_interface_id = near_interface_id_.DataField()[neighbour_index.first];
                            if (neighbor_near_interface_id[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]] <= -1)
                            {
                                Real phi_n_ = neighbor_phi[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                Vecd norm_to_face = neighbor_phi_gradient[neighbour_index.second[0]][neighbour_index.second[1]][neighbour_index.second[2]];
                                norm_to_face /= norm_to_face.norm() + TinyReal;
                                min_distance_n = SMIN(min_distance_n, (Vecd((Real)x, (Real)y, Real(z)) * data_spacing_ - phi_n_ * norm_to_face).norm());
                            }
                        });
                    phi_data[package_index][i][j][k] = min_distance_n;
                    // this immediate switch of near interface id
                    // does not intervening with the identification of unresolved interface
                    // based on the assumption that positive false_and negative bands are not close to each other
                    near_interface_id_data[package_index][i][j][k] = 1;
                }
            }
        });
}
//=================================================================================================//
void LevelSet::writeMeshFieldToPlt(std::ofstream &output_file)
{
    Arrayi number_of_operation = global_mesh_.AllGridPoints();

    output_file << "\n";
    output_file << "title='View'"
                << "\n";
    output_file << "variables= "
                << "x, "
                << "y, "
                << "z, "
                << "phi, "
                << "n_x, "
                << "n_y, "
                << "n_z "
                << "\n";
    output_file << "zone i=" << number_of_operation[0] << "  j=" << number_of_operation[1] << "  k=" << number_of_operation[2]
                << "  DATAPACKING=BLOCK  SOLUTIONTIME=" << 0 << "\n";

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[0] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[1] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                Vecd data_position = global_mesh_.GridPositionFromIndex(Arrayi(i, j, k));
                output_file << data_position[2] << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_, Arrayi(i, j, k))
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j, k))[0]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j, k))[1]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(phi_gradient_, Arrayi(i, j, k))[2]
                            << " ";
            }
            output_file << " \n";
        }

    for (int k = 0; k != number_of_operation[2]; ++k)
        for (int j = 0; j != number_of_operation[1]; ++j)
        {
            for (int i = 0; i != number_of_operation[0]; ++i)
            {
                output_file << DataValueFromGlobalIndex(near_interface_id_, Arrayi(i, j, k))
                            << " ";
            }
            output_file << " \n";
        }
}
//=============================================================================================//
RefinedLevelSet::RefinedLevelSet(BoundingBox tentative_bounds, LevelSet &coarse_level_set,
                                 Shape &shape, SPHAdaptation &sph_adaptation)
    : RefinedMesh(tentative_bounds, coarse_level_set, 4, shape, sph_adaptation)
{
    mesh_parallel_for(MeshRange(Arrayi::Zero(), all_cells_),
                      [&](size_t i, size_t j, size_t k)
                      {
                          initializeDataInACellFromCoarse(Arrayi(i, j, k));
                      });

    finishDataPackages();
}
//=============================================================================================//
} // namespace SPH
