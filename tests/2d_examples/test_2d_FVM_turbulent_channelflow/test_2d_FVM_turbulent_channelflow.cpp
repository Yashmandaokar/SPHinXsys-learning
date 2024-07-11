/**
 * @file 	2d_FVM_flow_around_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder in FVM.
 * @details We consider a flow passing by a cylinder in 2D in FVM framework.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#include "test_2d_FVM_turbulent_channelflow.h"
#include "common_weakly_compressible_FVM_classes.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh read_mesh_data(mesh_file_path);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, read_mesh_data.MinMeshEdge());
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, UnstructuredMesh>(ghost_boundary, read_mesh_data);
    
    GhostCreationFromMesh ghost_creation(water_block, read_mesh_data, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(water_block, read_mesh_data);
    SimpleDynamics<TCFInitialCondition> initial_condition(water_block);
    SimpleDynamics<fluid_dynamics::WallAdjacentCells> wall_adj_cell(water_block_inner, ghost_creation);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    //InteractionWithUpdate<fluid_dynamics::MeanVelocity> meanvelocity_relaxation(water_block_inner);
    //fluid_dynamics::WallAdjacentCells wall_adj_cell(water_block_inner, ghost_creation);
    InteractionWithUpdate<fluid_dynamics::KEpsilonStd1stHalf> tke(water_block_inner, ghost_creation);
    InteractionWithUpdate<fluid_dynamics::KEpsilonStd2ndHalf> dissipationrate(water_block_inner, ghost_creation);
    InteractionDynamics<fluid_dynamics::StdWallFunctionFVM> wall_function(water_block_inner, ghost_creation);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfInnerRiemannRANS> pressure_relaxation(water_block_inner, 0.0);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfInnerRiemann> density_relaxation(water_block_inner, 0.0);
    TCFBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(water_block, read_mesh_data.MinMeshEdge());
   // InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    //----------------------------------------------------------------------
//	Compute the force exerted on solid body due to fluid pressure and viscosity
//----------------------------------------------------------------------
    //InteractionDynamics<fluid_dynamics::ViscousForceFromFluidInFVM> viscous_force_on_solid(water_block_inner, ghost_creation.each_boundary_type_contact_real_index_);
    //InteractionDynamics<fluid_dynamics::PressureForceFromFluidInFVM<decltype(density_relaxation)>> pressure_force_on_solid(water_block_inner, ghost_creation.each_boundary_type_contact_real_index_);
   // visualization in FVM with data in cell
    BodyStatesRecordingInMeshToVtp write_real_body_states(water_block, read_mesh_data);
    //RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_on_inserted_body(water_block, "ViscousForceOnSolid");
    //ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_on_inserted_body(water_block, "PressureForceOnSolid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(water_block);

    initial_condition.exec();
    wall_adj_cell.exec();
    water_block_inner.updateConfiguration();
    water_block.addBodyStateForRecording<Vecd>("Velocity");
    water_block.addBodyStateForRecording<Real>("Density");
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("TKE");
    water_block.addBodyStateForRecording<Real>("Dissipation");
    water_block.addBodyStateForRecording<Real>("TurbulentViscosity");
    water_block.addBodyStateForRecording<Real>("TKEProduction");
    water_block.addBodyStateForRecording<Real>("TKEChangeRate");
    water_block.addBodyStateForRecording<Real>("DissipationChangeRate");
    water_block.addBodyStateForRecording<Real>("TKEAdvection");
    water_block.addBodyStateForRecording<Real>("TKELaplacian");
    water_block.addBodyStateForRecording<Real>("DissipationAdvection");
    water_block.addBodyStateForRecording<Real>("DissipationLaplacian");
    water_block.addBodyStateForRecording<Real>("DissipationProdscalar");
    water_block.addBodyStateForRecording<Real>("DissipationScalar");
    //water_block.addBodyStateForRecording<Vecd>("MeanVelocity");
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 20.0;
    Real output_interval = 0.1; /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            boundary_condition_setup.resetBoundaryConditions();
            tke.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            dissipationrate.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            wall_function.exec();
            boundary_condition_setup.resetBoundaryConditions();
            //viscous_force.exec();
            pressure_relaxation.exec(dt);
            write_maximum_speed.writeToFile(number_of_iterations);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
                    << GlobalStaticVariables::physical_time_
                    << "	dt = " << dt << "\n";
                //viscous_force_on_solid.exec();
                //pressure_force_on_solid.exec();
                //write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
                //write_total_pressure_force_on_inserted_body.writeToFile(number_of_iterations);
                //write_maximum_speed.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            if (number_of_iterations == 48)
            {
                Real k = 1.0;
            }

            write_real_body_states.writeToFile();
        }
        TickCount t2 = TickCount::now();
        //write_real_body_states.writeToFile();
        //write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;
    //write_total_viscous_force_on_inserted_body.testResult();
    return 0;
}
