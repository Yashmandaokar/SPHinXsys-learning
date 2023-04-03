/**
 * @file 	tank_v_multiphase.cpp
 * @brief 	3D Two-phase Sloshing in a Vertical Cylindrical Tank under Lateral Excitation
 */
#include "sphinxsys.h"
#include "tank_case.h"
/*@brief Namespace cite here.
*/
using namespace SPH;

/*
Main program starts here.
*/
int main(int ac, char* av[])
{
	/* Build up -- a SPHSystem -- */
	SPHSystem system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	//system.setRunParticleRelaxation(false);
	//// Tag for computation start with relaxed body fitted particles distribution.
	//system.setReloadParticles(true);
	/* Tag for computation from restart files. 0: start with initial condition. */
	system.setRestartStep(0);
	//handle command line arguments
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);

	/*
	@Brief creating body, materials and particles for the tank, water, and air.
	*/
	//SolidBody tank(system, makeShared<Tank>("Tank"));
	//tank.defineParticlesAndMaterial<SolidParticles, Solid>();
	///*tank.defineBodyLevelSetShape()->writeLevelSet(in_output);*/
	//(!system.RunParticleRelaxation() && system.ReloadParticles())
	//	? tank.generateParticles<ParticleGeneratorReload>(in_output, tank.getName())
	//	: tank.generateParticles<ParticleGeneratorLattice>();

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	/*water_block.defineBodyLevelSetShape()->writeLevelSet(in_output);*/
	water_block.generateParticles<ParticleGeneratorLattice>();


	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_a, c_f);
	/*air_block.defineBodyLevelSetShape()->writeLevelSet(in_output);*/
	air_block.generateParticles<ParticleGeneratorLattice>();
	

	//----------------------------------------------------------------------
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });
	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);

	
	SharedPtr<Gravity>gravity_ptr = makeShared<VariableGravity>();
	//SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block,gravity_ptr );
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, gravity_ptr);
	/* Fluid dynamics */
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> fluid_density_by_summation(water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_water_complex);
	
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> air_advection_time_step(air_block, U_g);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> air_acoustic_time_step(air_block);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemann> fluid_pressure_relaxation(water_air_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemann> fluid_density_relaxation(water_air_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemann> air_pressure_relaxation(air_water_complex);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemann> air_density_relaxation(air_water_complex);
		/** Confinement condition for wall and structure. */
	NearShapeSurface near_surface_water(water_block, makeShared<WallAndStructure>("WallAndStructure"));
	NearShapeSurface near_surface_air(air_block, makeShared<WallAndStructure>("WallAndStructure"));
	/*near_surface_water.level_set_shape_.writeLevelSet(in_output);*/
	fluid_dynamics::StaticConfinement confinement_condition_water(near_surface_water);
	fluid_dynamics::StaticConfinement confinement_condition_air(near_surface_air);
	fluid_density_by_summation.post_processes_.push_back(&confinement_condition_water.density_summation_);
	fluid_pressure_relaxation.post_processes_.push_back(&confinement_condition_water.pressure_relaxation_);
	fluid_density_relaxation.post_processes_.push_back(&confinement_condition_water.density_relaxation_);
	air_density_by_summation.post_processes_.push_back(&confinement_condition_air.density_summation_);
	air_pressure_relaxation.post_processes_.push_back(&confinement_condition_air.pressure_relaxation_);
	air_density_relaxation.post_processes_.push_back(&confinement_condition_air.density_relaxation_);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(in_output, system.real_bodies_);
	/*RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(in_output, water_block, gravity_ptr);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>;*/
		/*write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);*/
	/*BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("PorbeS1"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_1(in_output, probe_s1);
	BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_2(in_output, probe_s2);
	BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("PorbeS3"));
	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight, BodyRegionByCell>>
		probe_3(in_output, probe_s3);*/
	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	/*inner_normal_direction.parallel_exec();*/

	write_real_body_states.writeToFile(0);
	/*probe_1.writeToFile(0);
	probe_2.writeToFile(0);
	probe_3.writeToFile(0);*/
	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		/*water_block_contact.updateconfiguration();*/
		air_water_complex.updateConfiguration();
		/*air_block_contact.updateconfiguration();*/
	}
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real end_time = 5.0;		/**< End time. */
	Real output_interval = 0.1; /**< Time stamps for output of body states. */
	Real dt = 0.0;				/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile(0);
	/*write_water_mechanical_energy.writeToFile(0);*/
	/*write_recorded_water_pressure.writeToFile(0);*/
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_water_step.parallel_exec();
			Real Dt_f = fluid_advection_time_step.parallel_exec();
			Real Dt_a = air_advection_time_step.parallel_exec();
			Real Dt = SMIN(Dt_f, Dt_a);
			fluid_density_by_summation.parallel_exec();
			air_density_by_summation.parallel_exec();
			air_transport_correction.parallel_exec();
			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = fluid_acoustic_time_step.parallel_exec();
				Real dt_a = air_acoustic_time_step.parallel_exec();
				dt = SMIN(SMIN(dt_f, dt_a), Dt);
				fluid_pressure_relaxation.parallel_exec(dt);
				fluid_density_relaxation.parallel_exec(dt);
				air_pressure_relaxation.parallel_exec(dt);
				air_density_relaxation.parallel_exec(dt);
				/*dt = fluid_acoustic_time_step.parallel_exec();*/
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				//if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
				//{
				//	write_water_mechanical_energy.writeToFile(number_of_iterations);
				//	/*write_recorded_water_pressure.writeToFile(number_of_iterations);*/
				//}
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			water_block.updateCellLinkedListWithParticleSort(100);
			/*water_block_inner.updateConfiguration();*/
			air_block.updateCellLinkedListWithParticleSort(100);
			/*air_block_inner.updateConfiguration();*/
			air_water_complex.updateConfiguration();
			water_air_complex.updateConfiguration();
			/*fluid_observer_contact.updateConfiguration();*/
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}