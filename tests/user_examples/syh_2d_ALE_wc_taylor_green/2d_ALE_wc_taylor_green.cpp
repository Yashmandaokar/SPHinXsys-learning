/**
 * @file 	eulerian_taylor_green.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 2D eulerian_taylor_green vortex flow example.
 * @author 	Chi Zhang, Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
#include "ALE_weakly_compressible_fluid_dynamic_inner.h"
#include "ALE_weakly_compressible_fluid_dynamic_inner.hpp"

using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;					  /**< box length. */
Real DH = 1.0;					  /**< box height. */
Real resolution_ref = 1.0 / 100.0; /**< Global reference resolution. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					/**< Reference density of fluid. */
Real U_f = 1.0;						/**< Characteristic velocity. */
Real c_f = 10.0 * U_f;				/**< Reference sound speed. */
Real Re = 100;						/**< Reynolds number. */
Real mu_f = rho0_f * U_f * DL / Re; /**< Dynamics viscosity. */
Real heat_capacity_ratio = 1.4;		/**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> water_body_shape;
		water_body_shape.push_back(Vecd(0.0, 0.0));
		water_body_shape.push_back(Vecd(0.0, DH));
		water_body_shape.push_back(Vecd(DL, DH));
		water_body_shape.push_back(Vecd(DL, 0.0));
		water_body_shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(water_body_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Application dependent initial condition
//----------------------------------------------------------------------
class TaylorGreenInitialCondition
	: public eulerian_weakly_compressible_fluid_dynamics::WeaklyCompressibleFluidInitialCondition
{
public:
	explicit TaylorGreenInitialCondition(SPHBody &sph_body)
		: eulerian_weakly_compressible_fluid_dynamics::WeaklyCompressibleFluidInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		/** initial momentum and energy profile */
		rho_[index_i] = rho0_f;
		/*p_[index_i] = pow(c_f, 2) * rho_[index_i] / gamma_;*/
		vel_[index_i][0] = -cos(2.0 * Pi * pos_[index_i][0]) *
						   sin(2.0 * Pi * pos_[index_i][1]);
		vel_[index_i][1] = sin(2.0 * Pi * pos_[index_i][0]) *
						   cos(2.0 * Pi * pos_[index_i][1]);
		mom_[index_i] = rho_[index_i] * vel_[index_i];
		/*Real rho_e = p_[index_i] / (gamma_ - 1.0);*/
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	EulerianFluidBody water_body(sph_system, makeShared<WaterBlock>("WaterBody"));
	water_body.defineParticlesAndMaterial<WeaklyCompressibleFluidParticles, WeaklyCompressibleFluid>(rho0_f, heat_capacity_ratio, mu_f);
	water_body.generateParticles<ParticleGeneratorLattice>();
	water_body.addBodyStateForRecording<Real>("Density");
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_body_inner(water_body);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initial condition with momentum and energy field */
	SimpleDynamics<TaylorGreenInitialCondition> initial_condition(water_body);
	/** Initialize a time step. */
	SimpleDynamics<ALE_weakly_compressible_fluid_dynamic::ALEFlowTimeStepInitialization> time_step_initialization(water_body);
	/** Computing viscous acceleration. */
	Dynamics1Level<ALE_weakly_compressible_fluid_dynamic::Integration1stHalfHLLCWithLimiterRiemann> pressure_relaxation(water_body_inner);
	InteractionWithUpdate<ALE_weakly_compressible_fluid_dynamic::Integration2ndHalfHLLCWithLimiterRiemann> density_relaxation(water_body_inner);
	InteractionWithUpdate<ALE_weakly_compressible_fluid_dynamic::DensitySummationInner> update_density_by_summation(water_body_inner);
	InteractionDynamics<ALE_weakly_compressible_fluid_dynamic::ViscousAccelerationInner> viscous_acceleration(water_body_inner);
	InteractionDynamics<ALE_weakly_compressible_fluid_dynamic::TransportVelocityCorrectionInner> transport_velocity_correction(water_body_inner);
	ReduceDynamics<ALE_weakly_compressible_fluid_dynamic::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_body, U_f);
	ReduceDynamics<ALE_weakly_compressible_fluid_dynamic::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	PeriodicConditionUsingCellLinkedList periodic_condition_x(water_body, water_body.getBodyShapeBounds(), xAxis);
	PeriodicConditionUsingCellLinkedList periodic_condition_y(water_body, water_body.getBodyShapeBounds(), yAxis);
	/** Time step size with considering sound wave speed. */
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	ReloadParticleIO write_particle_reload_files(io_environment, water_body);
	RegressionTestEnsembleAveraged<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_total_mechanical_energy(io_environment, water_body);
	/** Output the maximum speed of the fluid body. */
	RegressionTestEnsembleAveraged<ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>>>
		write_maximum_speed(io_environment, water_body);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	initial_condition.exec();
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_x.update_cell_linked_list_.exec();
	periodic_condition_y.update_cell_linked_list_.exec();
	sph_system.initializeSystemConfigurations();
	
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real end_time = 5.0;
	Real output_interval = 0.1; /**< Time stamps for output of body states. */
	Real dt = 0;
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	body_states_recording.writeToFile();
	/** Output the mechanical energy of fluid. */
	write_total_mechanical_energy.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			time_step_initialization.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();
			update_density_by_summation.exec();
			viscous_acceleration.exec();
			transport_velocity_correction.exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				// avoid possible smaller acoustic time step size for viscous flow
				dt = SMIN(get_fluid_time_step_size.exec(), Dt);
				relaxation_time += dt;
				integration_time += dt;
				pressure_relaxation.exec(dt);
				density_relaxation.exec(dt);
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			/** Water block configuration and periodic condition. */
			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			water_body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			water_body_inner.updateConfiguration();
		}

		TickCount t2 = TickCount::now();
		write_total_mechanical_energy.writeToFile(number_of_iterations);
		write_maximum_speed.writeToFile(number_of_iterations);
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}

	    TickCount t2 = TickCount::now();
		
		/*write_total_mechanical_energy.writeToFile(number_of_iterations);
		write_maximum_speed.writeToFile(number_of_iterations);*/
	
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	
		TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		 << " seconds." << std::endl;

	/*write_total_mechanical_energy.newResultTest();
	write_maximum_speed.newResultTest();*/

	return 0;
}