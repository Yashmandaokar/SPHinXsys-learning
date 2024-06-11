/**
 * @file 	2d_eulerian_flow_around_cylinder_LG.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder coupling with the Laguerre Gauss kernel.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#ifndef TEST_2D_EULERIAN_TURBULENT_CHANNELFLOW_H
#define TEST_2D_EULERIAN_TURBULENT_CHANNELFLOW_H
#include "sphinxsys.h"
#include "TurbulenceModel.h"
#include "rans_fluid_integration.hpp"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 1;                                           /**< Channel height. */
//Real DL = 120 * DH / 2.0;                             /**< Channel length. */
Real DL = 8 * M_PI * DH / 2.0;
Real resolution_ref = DH / 6.0;                     /**< Initial reference particle spacing. */
Real BW = resolution_ref;                       /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;
Real U_f = 1.0;
Real rey_bulk = 5186;
Real mu_f = (rho0_f * U_f * DH * 0.5) / rey_bulk;       /**< Dynamic Viscosity. */
Real c_f = 10.0 * U_f;                               /**< Reference sound speed. */
BoundingBox system_domain_bounds(Vecd(0, -BW), Vecd(DL, DH + BW));
Real gravity_g = 12.0 * mu_f * U_f / rho0_f / DH / DH; /**< Gravity force of fluid. */

Real I = 0.05;
Real viscosityRatio = 2;
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------

std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));

    return water_block_shape;
}
Vecd channel_halfsize = Vecd(0.5 * DL, 0.5 * DH);
Vecd channel_translation = Vecd(-DL, 0.0) + channel_halfsize;

class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
    }
};
/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(0.0, -BW));
        outer_wall_shape.push_back(Vecd(0.0, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, -BW));
        outer_wall_shape.push_back(Vecd(0.0, -BW));
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(0.0, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, DH));
        inner_wall_shape.push_back(Vecd(DL, DH));
        inner_wall_shape.push_back(Vecd(DL, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class TCFInitialCondition
    : public fluid_dynamics::FluidInitialCondition
{
public:
    explicit TCFInitialCondition(SPHBody& sph_body)
        : FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
        rho_(particles_->rho_), mass_(particles_->mass_), Vol_(particles_->Vol_),
        p_(*particles_->getVariableByName<Real>("Pressure")), C_mu_(0.09)
    {
        particles_->registerVariable(mom_, "Momentum");
        particles_->registerVariable(K_, "TKE");
        particles_->registerVariable(Eps_, "Dissipation");
        particles_->registerVariable(mu_t_, "TurbulentViscosity");
    };

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = rho0_f;
        p_[index_i] = 50.0 / 117.6655;
        vel_[index_i][0] = 0.5;
        vel_[index_i][1] = 0.0;
        K_[index_i] = (3.0 / 2.0) * (vel_[index_i].squaredNorm()) * (I * I);
        Eps_[index_i] = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (mu_f)) * (1 / viscosityRatio);
    }
protected:
    StdLargeVec<Vecd>& pos_, &vel_;
    StdLargeVec<Real>& rho_, & mass_, & Vol_, & p_;
    StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
    StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
    Real C_mu_;
    StdLargeVec<Real> mu_t_, K_, Eps_;
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape& aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType& boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd& position, Vecd& velocity)
    {
        Vecd target_velocity = velocity;
        target_velocity[0] = 1.0;
        return target_velocity;
    }
};
#endif //TEST_2D_EULERIAN_TURBULENT_CHANNELFLOW_H
