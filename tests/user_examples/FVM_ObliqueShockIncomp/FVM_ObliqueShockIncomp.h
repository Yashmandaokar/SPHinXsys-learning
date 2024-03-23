/**
 * @file 	Yash_3d_FVM_double_mach_reflection_Unstruct.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_3D_Unstruct_H
#define FVM_3D_Unstruct_H
#include "unstructured_mesh.h"              // shared eulerian classes for weakly-compressible and compressible fluid in FVM.
#include "common_weakly_compressible_FVM_classes_3D.hpp" // classes for weakly compressible fluid only in FVM.
#include "eulerian_fluid_dynamics.hpp"                // eulerian classes for weakly compressible fluid only.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0; //0.769846;                            /**< Computation domain length. */
Real DH = 0.6494805454;//0.5;                           /**< Computation domain height. */
Real DW = 0.038968832;//0.03;                          /**< Computation domain width. */
Real particle_spacing_ref = 1.0 / 500.0; /**< Initial reference particle spacing. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-0.38968832, 0.0, 0.0), Vec3d(DL, DH, DW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
//Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = 0; //rho0_f * U_f * (DL) / Re;              /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_fullpath = "./input/Channel.msh";

//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------

class WaveBody : public ComplexShape
{
public:
    explicit WaveBody(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wave(0.5 * DH, 0.5 * DL, 0.5 * DW);
        Transform translation_wave(halfsize_wave);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wave), halfsize_wave);
    }
};
//______________________________________________________________________
// Initialization
//______________________________________________________________________

///----------------------------------------------------------------------
//	DMFBoundaryConditionSetup
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public fluid_dynamics::FluidDataInner
{
public:
    FACBoundaryConditionSetup(BaseInnerRelationInFVM& inner_relation, vector<vector<size_t>> each_boundary_type_with_all_ghosts_index,
        vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_, vector<vector<size_t>> each_boundary_type_contact_real_index)
        : fluid_dynamics::FluidDataInner(inner_relation), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
        Vol_(particles_->Vol_), vel_(particles_->vel_), mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
        fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())), total_ghost_particles_(particles_->total_ghost_particles_),
        real_particles_bound_(particles_->real_particles_bound_), each_boundary_type_with_all_ghosts_index_(each_boundary_type_with_all_ghosts_index),
        each_boundary_type_with_all_ghosts_eij_(each_boundary_type_with_all_ghosts_eij_), each_boundary_type_contact_real_index_(each_boundary_type_contact_real_index) {};
    virtual ~FACBoundaryConditionSetup() {};

    void resetBoundaryConditions()
    {
        for (size_t boundary_type = 0; boundary_type < each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
        {
            if (!each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
            {
                for (size_t ghost_number = 0; ghost_number != each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                {
                    size_t ghost_index = each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                    size_t index_i = each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                    Vecd e_ij = each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];
                    if (boundary_type == 3)
                    {
                        // non-slip wall boundary
                        vel_[ghost_index] = (vel_[index_i] - e_ij.dot(vel_[index_i]) * (e_ij)) + (-e_ij.dot(vel_[index_i]) * (e_ij));
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                    }
                    // Velocity Inlet boundary
                    if (boundary_type == 10)
                    {
                        Vecd far_field_velocity(1.0, 0.0, 0.0);
                        vel_[index_i] = far_field_velocity;
                        vel_[ghost_index] = far_field_velocity;
                        p_[ghost_index] = p_[index_i];
                        rho_[ghost_index] = rho_[index_i];
                    }
                    if (boundary_type == 5)
                    {
                        // Pressure Outlet boundary condition
                        vel_[ghost_index] = vel_[index_i];
                        p_[index_i] = 1.0;
                        p_[ghost_index] = 1.0;
                        rho_[index_i] = rho0_f;
                        rho_[ghost_index] = rho_[index_i];
                    }

                    // SYMMETRY boundary condition
                    if (boundary_type == 7)
                    {
                        vel_[ghost_index] = (vel_[index_i] - 2 * e_ij.dot(vel_[index_i]) * e_ij);
                        rho_[ghost_index] = rho_[index_i];
                        p_[ghost_index] = p_[index_i];

                    }
                }
            }
        }
    };

protected:
    StdLargeVec<Real>& rho_, & p_, & Vol_;
    StdLargeVec<Vecd>& vel_, & mom_, & pos_;
    Fluid& fluid_;
    size_t& total_ghost_particles_;
    size_t& real_particles_bound_;
    vector<vector<size_t>> each_boundary_type_with_all_ghosts_index_;
    vector<vector<Vecd>> each_boundary_type_with_all_ghosts_eij_;
    vector<vector<size_t>> each_boundary_type_contact_real_index_;
};
#endif // FVM_3D_Unstruct_H