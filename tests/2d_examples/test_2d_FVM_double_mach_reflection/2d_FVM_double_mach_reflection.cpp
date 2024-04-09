/**
 * @file 	2d_FVM_double_mach_reflection.cpp
 * @brief 	This is the compressible test for the realization of FVM in the SPHinXsys.
 * @details We consider a double mach reflection case.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "2d_FVM_double_mach_reflection.h"
#include "sphinxsys.h"

using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh ansys_mesh(double_mach_reflection_mesh1_fullpath);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, ansys_mesh.MinMeshEdge());
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody wave_block(sph_system, makeShared<WaveBody>("WaveBody"));
    wave_block.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_another, heat_capacity_ratio);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    wave_block.generateParticlesWithReserve<UnstructuredMesh>(ghost_boundary, ansys_mesh);
    wave_block.addBodyStateForRecording<Real>("Density");
    wave_block.addBodyStateForRecording<Real>("Pressure");
    /** Initial condition and register variables*/
    SimpleDynamics<DMFInitialCondition> initial_condition(wave_block);
    GhostCreationFromMesh ghost_creation(wave_block, ansys_mesh, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(wave_block, ansys_mesh);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Boundary conditions set up */
    DMFBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<CompressibleAcousticTimeStepSizeInFVM> get_fluid_time_step_size(wave_block, ansys_mesh.MinMeshEdge(), 0.2);
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfHLLCRiemann> pressure_relaxation(water_block_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfHLLCRiemann> density_relaxation(water_block_inner);
    // Visualization in FVM with date in cell.
    BodyStatesRecordingInMeshToVtp write_real_body_states(wave_block, ansys_mesh);
    RegressionTestEnsembleAverage<ReducedQuantityRecording<MaximumSpeed>>
        write_maximum_speed(wave_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_block_inner.updateConfiguration();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 0.2;
    Real output_interval = 0.01; /**< time stamps for output. */
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
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);

            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                cout << fixed << setprecision(9) << "N=" << number_of_iterations << "	Time = "
                     << GlobalStaticVariables::physical_time_
                     << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

    if (sph_system.GenerateRegressionData())
    {
        write_maximum_speed.generateDataBase(1.0e-3, 1.0e-3);
    }
    else
    {
        write_maximum_speed.testResult();
    }
    return 0;
}
