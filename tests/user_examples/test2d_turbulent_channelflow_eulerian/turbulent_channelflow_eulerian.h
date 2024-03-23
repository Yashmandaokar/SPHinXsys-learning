/**
 * @file 	2d_eulerian_flow_around_cylinder_LG.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder coupling with the Laguerre Gauss kernel.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#ifndef TURBULENT_CHANNELFLOW_EULERIAN_H
#define TURBULENT_CHANNELFLOW_EULERIAN_H
#include "eulerian_fluid_dynamics.hpp" // eulerian classes for weakly compressible fluid only.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2;                        /**< Channel length. */
Real DH = 0.5;                        /**< Channel height. */
Real resolution_ref = DH / 20.0;       /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;    /**< Extending width for BCs. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                        /**< Density. */
Real Re = 100.0;
Real gravity_g = 1.0e-4;                                 /**< Gravity force of fluid. */
Real mu_f = 1.0e-6;                                     /**< Viscosity. */
Real U_f = gravity_g * DH * DH / mu_f;                  /**< Characteristic velocity. */           
Real c_f = 10.0 * U_f;                                  /**< Speed of sound. */     
//Real mu_f = rho0_f * U_f * (2.0 * DH) / Re;           /**< Dynamics viscosity. */

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
        outer_wall_shape.push_back(Vecd(0, -BW));
        outer_wall_shape.push_back(Vecd(0, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, DH + BW));
        outer_wall_shape.push_back(Vecd(DL, -BW));
        outer_wall_shape.push_back(Vecd(0, -BW));
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
class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
public:
    explicit FarFieldBoundary(BaseInnerRelation& inner_relation)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation)
    {
        rho_farfield_ = rho0_f;
        sound_speed_ = c_f;
        vel_farfield_ = Vecd(U_f, 0.0);
    };
    virtual ~FarFieldBoundary() {};
};
#endif //TURBULENT_CHANNELFLOW_EULERIAN_H
