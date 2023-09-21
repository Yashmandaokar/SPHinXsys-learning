/**
 * @file 	tank_v_multiphase.cpp
 * @brief 	3D Two-phase Sloshing in a Vertical Cylindrical Tank under Lateral Excitation
 */
#include "sphinxsys.h"
#include "YM_tank_case.h"
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
	system.handleCommandlineOptions(ac, av);
	/* Output environment. */
	IOEnvironment in_output(system);

	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.setRunParticleRelaxation(true);
	/** Tag for starting with relaxed body-fitted particles distribution */
	system.setReloadParticles(false);
	/*
	@Brief creating body, materials and particles for the tank, water, and air.
	*/

	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
	//(!system.RunParticleRelaxation() && system.ReloadParticles())
		//? water_block.generateParticles<ParticleGeneratorReload>(in_output, water_block.getName())
	water_block.generateParticles<ParticleGeneratorLattice>();
	water_block.addBodyStateForRecording<Vecd>("Acceleration");
	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<Real>("Density");


	FluidBody air_block(system, makeShared<AirBlock>("AirBody"));
	air_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_a, c_f);
	air_block.defineBodyLevelSetShape()->writeLevelSet(in_output);
//	(!system.RunParticleRelaxation() && system.ReloadParticles())
//		? air_block.generateParticles<ParticleGeneratorReload>(in_output, air_block.getName())
	air_block.generateParticles<ParticleGeneratorLattice>();
	air_block.addBodyStateForRecording<Real>("Pressure");

	SolidBody tank(system, makeShared<Tank>("Tank"));
	tank.defineParticlesAndMaterial<SolidParticles, Solid>();
	tank.defineBodyLevelSetShape()->writeLevelSet(in_output);
//	tank.defineComponentLevelSetShape("OuterWall");
//	tank.defineComponentLevelSetShape("InnerWall");
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? tank.generateParticles<ParticleGeneratorReload>(in_output, tank.getName())
		: tank.generateParticles<ParticleGeneratorLattice>();
	tank.addBodyStateForRecording<Vecd>("NormalDirection");

//	ObserverBody fluid_observer(system, "FluidObserver");
//	fluid_observer.generateParticles<WaterObserverParticleGenerator>();

	//--------------------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//--------------------------------------------------------------------------------

	InnerRelation tank_inner(tank);
	InnerRelation water_block_inner(water_block);
	InnerRelation air_block_inner(air_block);

	ComplexRelation tank_boundary_complex(tank, RealBodyVector{ &water_block, &air_block });
	ComplexRelation water_air_complex(water_block, { &air_block });
	ComplexRelation air_water_complex(air_block, { &water_block });

	ContactRelation water_tank_contact(water_block, { &tank });
	ContactRelation air_tank_contact(air_block, { &tank });
	ContactRelation tank_boundary_water_contact(tank, RealBodyVector{ &water_block,&air_block });



	//	ContactRelation fluid_observer_contact(fluid_observer, { &water_block });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the tank particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_tank_particles(tank);
//		SimpleDynamics<RandomizeParticlePosition> random_water_particles(water_block);
//		SimpleDynamics<RandomizeParticlePosition> random_air_particles(air_block);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_tank_to_vtp(in_output, { &tank });
//		BodyStatesRecordingToVtp write_water_to_vtp(in_output, { &water_block });
//		BodyStatesRecordingToVtp write_air_to_vtp(in_output, { &air_block });
		/** Write the particle reload files. */
		ReloadParticleIO write_tank_boundary_particle_reload_files(in_output, tank, "Tank");
//		ReloadParticleIO write_water_particle_reload_files(in_output, water_block, "WaterBody");
//		ReloadParticleIO write_air_particle_reload_files(in_output, air_block, "AirBody");
		/** A  Physics relaxation step. */
//		relax_dynamics::RelaxationStepInner water_relaxation_step_inner(water_block_inner);
//		relax_dynamics::RelaxationStepInner air_relaxation_step_inner(air_block_inner);
		relax_dynamics::RelaxationStepInner tank_relaxation_step_inner(tank_inner);
//		relax_dynamics::RelaxationStepComplex tank_boundary_relaxation_step_complex(tank_boundary_complex, "Tank",true);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_tank_particles.exec(0.25);
	//	random_water_particles.exec(0.25);
	//	random_air_particles.exec(0.25);
		tank_relaxation_step_inner.SurfaceBounding().exec();
	//	tank_boundary_relaxation_step_complex.SurfaceBounding().exec();
	//	water_relaxation_step_inner.SurfaceBounding().exec();
	//	air_relaxation_step_inner.SurfaceBounding().exec();
		write_tank_to_vtp.writeToFile(0);
	//	write_water_to_vtp.writeToFile(0);
	//	write_air_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the tank.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			tank_relaxation_step_inner.exec();
		//	water_relaxation_step_inner.exec();
		//	air_relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the tank N = " << ite_p << "\n";
			//	write_water_to_vtp.writeToFile(ite_p);
			//	write_air_to_vtp.writeToFile(ite_p);
				write_tank_to_vtp.writeToFile(ite_p);

			}
		}
		std::cout << "The physics relaxation process of finish !" << std::endl;
		/** Output results. */
		write_tank_boundary_particle_reload_files.writeToFile(0);
	//	write_water_particle_reload_files.writeToFile(0);
	//	write_air_particle_reload_files.writeToFile(0);

		return 0;
	}

	/*
	@Brief define simple data file input and outputs functions.
	*/
	BodyStatesRecordingToVtp 			write_real_body_states(in_output, system.real_bodies_);
	RestartIO							restart_io(in_output, system.real_bodies_);
	// --------------------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//--------------------------------------------------------------------------------
	InteractionWithUpdate<CorrectedConfigurationInner> tank_corrected_configuration(tank_inner);
	SimpleDynamics<NormalDirectionFromShapeAndOp> inner_normal_direction(tank,"InnerWall");
//	InteractionDynamics<solid_dynamics::CorrectConfiguration> tank_boundary_corrected_configuration(tank_boundary_complex.getInnerRelation());
	SimpleDynamics<TimeStepInitialization> initialize_a_water_step(water_block, makeShared<VariableGravity>());
	SimpleDynamics<TimeStepInitialization> initialize_a_air_step(air_block, makeShared<VariableGravity>());

	//SharedPtr<Gravity>gravity_ptr = makeShared<VariableGravity>();

	//--------------------------------------------------------------------------------
	//	Algorithms of fluid dynamics.
	//--------------------------------------------------------------------------------
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> 
		update_water_density_by_summation(water_tank_contact, water_air_complex.getInnerRelation());
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_air_density_by_summation(air_tank_contact, air_water_complex);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_tank_contact, air_water_complex);

	/*InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> air_density_by_summation(air_block_contact,air_water_complex);*/
	/*InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex> air_transport_correction(air_block_contact,air_water_complex);*/
	/** Time step size without considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_water_advection_time_step(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_air_advection_time_step(air_block, U_g);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_water_acoustic_time_step(water_block);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_air_acoustic_time_step(air_block);

	/** Riemann slover for pressure and density relaxation */
//	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfRiemannWithWall> water_pressure_relaxation(water_tank_contact, water_air_complex);
//	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall> water_density_relaxation(water_tank_contact, water_air_complex);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> water_pressure_relaxation(water_tank_contact, water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemannWithWall> water_density_relaxation(water_tank_contact, water_air_complex.getInnerRelation());
	Dynamics1Level<fluid_dynamics::ExtendMultiPhaseIntegration1stHalfRiemannWithWall> air_pressure_relaxation(air_tank_contact, air_water_complex, 2);
	Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfRiemannWithWall> air_density_relaxation(air_tank_contact, air_water_complex);

//	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_water(water_tank_contact, water_air_complex);
//	InteractionDynamics<fluid_dynamics::ViscousAccelerationMultiPhaseWithWall> viscous_acceleration_air(air_tank_contact, air_water_complex);
	
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	/*RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(in_output, water_block, gravity_ptr);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>;
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);*/
//	BodyRegionByCell probe_s1(water_block, makeShared<ProbeS1>("PorbeS1"));
//	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
//		probe_1(in_output, probe_s1);
//	BodyRegionByCell probe_s2(water_block, makeShared<ProbeS2>("PorbeS2"));
//	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
//		probe_2(in_output, probe_s2);
//	BodyRegionByCell probe_s3(water_block, makeShared<ProbeS3>("PorbeS3"));
//	ReducedQuantityRecording<ReduceDynamics<fluid_dynamics::FreeSurfaceHeight>>
//		probe_3(in_output, probe_s3);
	/**
 * @brief Pre-simulation.
 */
 /** initialize cell linked lists for all bodies. */
	system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	system.initializeSystemConfigurations();
	/** computing surface normal direction for the tank. */
	tank_corrected_configuration.exec();
	inner_normal_direction.exec();


	write_real_body_states.writeToFile(0);

	if (system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.RestartStep());
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		water_tank_contact.updateConfiguration();				
		air_water_complex.updateConfiguration();
		air_tank_contact.updateConfiguration();
		//fluid_observer_contact.updateConfiguration();
	}
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 16.0;		/**< End time. */
	Real output_interval = 0.025; /**< Time stamps for output of body states. */
	Real dt = 0.0;				/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	TickCount::interval_t interval_computing_time_step;
	TickCount::interval_t interval_computing_pressure_relaxation;
	TickCount::interval_t interval_updating_configuration;
	TickCount time_instance;
	
	//write_water_mechanical_energy.writeToFile(0);
	/*write_recorded_water_pressure.writeToFile(0);*/
	write_real_body_states.writeToFile(0);
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
			time_instance = TickCount::now();
			initialize_a_water_step.exec();
			initialize_a_air_step.exec();

			Real Dt_f = get_water_advection_time_step.exec();
			Real Dt_a = get_air_advection_time_step.exec();
			Real Dt = SMIN(Dt_f, Dt_a);

			update_water_density_by_summation.exec();
			update_air_density_by_summation.exec();
			air_transport_correction.exec();

//			viscous_acceleration_air.exec();
//			viscous_acceleration_water.exec();

			interval_computing_time_step += TickCount::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = TickCount::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_acoustic_time_step.exec();
				Real dt_a = get_air_acoustic_time_step.exec();

				dt = SMIN(SMIN(dt_f, dt_a), Dt);
				water_pressure_relaxation.exec(dt);
				air_pressure_relaxation.exec(dt);

				water_density_relaxation.exec(dt);
				air_density_relaxation.exec(dt);


				size_t inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += TickCount::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";


				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = TickCount::now();

			water_block.updateCellLinkedListWithParticleSort(100);
			air_block.updateCellLinkedListWithParticleSort(100);
			tank.updateCellLinkedList();
			air_water_complex.updateConfiguration();
			air_tank_contact.updateConfiguration();
			water_air_complex.updateConfiguration();
			water_tank_contact.updateConfiguration();
//			tank_boundary_water_contact.updateConfiguration();
			//fluid_observer_contact.updateConfiguration();
			interval_updating_configuration += TickCount::now() - time_instance;
			//if (number_of_iterations != 0 && number_of_iterations % observation_sample_interval == 0)
			//{
				//write_water_mechanical_energy.writeToFile(number_of_iterations);
				/*write_recorded_water_pressure.writeToFile(number_of_iterations);*/
			//}
		}

		TickCount t2 = TickCount::now();
		write_real_body_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;

	}
	TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
		<< interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";
	return 0;
}
