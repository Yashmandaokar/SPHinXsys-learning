/**
 * @file 	Yash_3d_FVM_double_mach_reflection_Unstruct.h
 * @brief 	This is a test to show the double mach reflection in FVM.
 * @details See https://doi.org/10.1016/j.jcp.2010.08.019 for the detailed problem setup.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef FVM_3D_Unstruct_OLD_H
#define FVM_3D_Unstruct_OLD_H
#include "unstructured_mesh.h"              // shared eulerian classes for weakly-compressible and compressible fluid in FVM.
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                            /**< Computation domain length. */
Real DH = 1.0;                           /**< Computation domain height. */
Real DW = 0.1;                          /**< Computation domain width. */
Real particle_spacing_ref = 1.0 / 240.0; /**< Initial reference particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(0.0, 0.0, 0.0), Vec3d(DL, DH, DW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_one = 1.4;                         /**< initial density of one fluid. */
Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
Real w_one = 0.0;							/**< initial velocity of one fluid in Z axis. */
Real p_one = 1.0;                            /**< initial pressure of one fluid. */
Real rho0_another = 8.0;                     /**< initial density of another. */
Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
Real w_another = 1; /**< initial velocity of another in Z axis. */
Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string double_mach_reflection_Unstruct_mesh_fullpath = "./input/double_mach_reflection_3D_fluent_ref.msh";
//
//	Define geometries and body shapes
//----------------------------------------------------------------------
/*
std::vector<Vecd> CreatComputationDomain()
{
    // geometry
    std::vector<Vecd> computation_domain;
    computation_domain.push_back(Vecd(0.0, 0.0, 0.0));
    computation_domain.push_back(Vecd(0.0, DH, 0.0));
    computation_domain.push_back(Vecd(DL, DH, 0.0));
    computation_domain.push_back(Vecd(DL, 0.0, 0.0));
    computation_domain.push_back(Vecd(0.0, 0.0, 0.0));
    computation_domain.push_back(Vecd(0.0, 0.0, DW));
    computation_domain.push_back(Vecd(0.0, DH, DW));
    computation_domain.push_back(Vecd(DL, DH, DW));
    computation_domain.push_back(Vecd(DL, 0.0, DW));
    computation_domain.push_back(Vecd(0.0, 0.0, DW));

    return computation_domain;
}*/

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
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DMFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit DMFInitialCondition(SPHBody& sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
        rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
        particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
        particles_->registerVariable(E_, "TotalEnergy");
        particles_->registerVariable(dE_dt_, "TotalEnergyChangeRate");
        particles_->registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
        gamma_ = heat_capacity_ratio;
    };
    void update(size_t index_i, Real dt)
    {
        if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - 1.0 / 6.0))
        {
            /** initial left wave pressure,momentum and energy profile */
            rho_[index_i] = rho0_another;
            p_[index_i] = p_another;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_another;
            vel_[index_i][1] = v_another;
            vel_[index_i][2] = w_another;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
        else
        {
            rho_[index_i] = rho0_one;
            p_[index_i] = p_one;
            Real rho_e = p_[index_i] / (gamma_ - 1.0);
            vel_[index_i][0] = u_one;
            vel_[index_i][1] = v_one;
            vel_[index_i][2] = w_another;
            mom_[index_i] = rho_[index_i] * vel_[index_i];
            E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
        }
    }

protected:
    StdLargeVec<Vecd>& pos_, & vel_;
    StdLargeVec<Real>& rho_, & p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real gamma_;
};
#endif // FVM_3D_Unstruct_OLD_H