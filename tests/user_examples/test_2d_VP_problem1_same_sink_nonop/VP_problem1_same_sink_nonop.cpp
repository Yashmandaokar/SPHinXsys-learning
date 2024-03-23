/**
 * @file 	VP_problem1_non_optimized
 * @brief 	This is the steady test for the same sink (2/10) .
 * @author 	Bo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
#include <gtest/gtest.h>
using namespace SPH; //Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 50.0;
Real BW = resolution_ref * 4.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1;
std::array<std::string, 1> species_name_list{ "Phi" };
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_ = 0.0;
Real high_ = 300.0;
Real low_ = 300.0;
Real heat_source = 1000.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
	std::vector<Vecd> thermalDomainShape;
	thermalDomainShape.push_back(Vecd(0.0, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, H));
	thermalDomainShape.push_back(Vecd(L, H));
	thermalDomainShape.push_back(Vecd(L, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, 0.0));

	return thermalDomainShape;
};

std::vector<Vecd> createBoundaryDomain()
{
	std::vector<Vecd> boundaryDomain;
	boundaryDomain.push_back(Vecd(-BW, -BW));
	boundaryDomain.push_back(Vecd(-BW, H + BW));
	boundaryDomain.push_back(Vecd(L + BW, H + BW));
	boundaryDomain.push_back(Vecd(L + BW, -BW));
	boundaryDomain.push_back(Vecd(-BW, -BW));

	return boundaryDomain;
};

//----------------------------------------------------------------------
//	Define SPH bodies. 
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
public:
	explicit DiffusionBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
	}
};

class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createBoundaryDomain(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::sub);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties. 
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
public:
	DiffusionMaterial() : DiffusionReaction<Solid>({ "Phi" }, SharedPtr<NoReaction>())
	{
		initializeAnDiffusion<LocalIsotropicDiffusion>("Phi", "Phi", diffusion_coeff);
	}
};

using DiffusionParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
using WallParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
//----------------------------------------------------------------------
//	Application dependent initial condition. 
//----------------------------------------------------------------------
class DiffusionBodyInitialCondition
	: public DiffusionReactionInitialCondition<DiffusionParticles>
{
protected:
	size_t phi_;
	StdLargeVec<Real>& heat_source_;

public:
	explicit DiffusionBodyInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialCondition<DiffusionParticles>(sph_body),
		heat_source_(*(particles_->getVariableByName<Real>("HeatSource")))
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = 550 + 50 * (double)rand() / RAND_MAX;
		heat_source_[index_i] = heat_source;
	};
};

class WallBoundaryInitialCondition
	: public DiffusionReactionInitialCondition<WallParticles>
{
protected:
	size_t phi_;

public:
	explicit WallBoundaryInitialCondition(SolidBody& diffusion_body)
		: DiffusionReactionInitialCondition<WallParticles>(diffusion_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	}

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = -0.0;
		if (pos_[index_i][1] < 0 && pos_[index_i][00] > 0.4 * L && pos_[index_i][0] < 0.6 * L)
		{
			all_species_[phi_][index_i] = low_;
		}
		if (pos_[index_i][1] > 1 && pos_[index_i][0] > 0.4 * L && pos_[index_i][0] < 0.6 * L)
		{
			all_species_[phi_][index_i] = high_;
		}
	};
};
//----------------------------------------------------------------------
//	An observer body to measure  at given positions. 
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	ObserverParticleGenerator(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		/** A line of measuring points at the middle line. */
		size_t number_of_observation_points = 10;
		Real range_of_measure = L;
		Real start_of_measure = 0;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(0.5 * L, range_of_measure * (Real)i /
				(Real)(number_of_observation_points - 1) + start_of_measure);
			positions_.push_back(point_coordinate);
		}
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
TEST(test_optimization, test_problem1_non_optimized)
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
	diffusion_body.defineParticlesAndMaterial<DiffusionParticles, DiffusionMaterial>();
	diffusion_body.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<WallParticles, DiffusionMaterial>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	//----------------------------  ------------------------------------------
	//	Particle and body creation of  observers.
	//----------------------------------------------------------------------
	ObserverBody _observer(sph_system, "Observer");
	_observer.generateParticles<ObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the range of bodies to build neighbor particle lists.
	ComplexRelation diffusion_body_complex(diffusion_body, { &wall_boundary });
	ContactRelation _observer_contact(_observer, { &diffusion_body });
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	InteractionSplit<SplittingByPDEWithBoundary<DiffusionParticles, WallParticles, Real>>
		_splitting(diffusion_body_complex, "Phi");
	GetDiffusionTimeStepSize<DiffusionParticles> get_time_step_size(diffusion_body);
	SimpleDynamics<DiffusionBodyInitialCondition> setup_diffusion_initial_condition(diffusion_body);
	SimpleDynamics<WallBoundaryInitialCondition> setup_boundary_condition(wall_boundary);
	ReduceAverage<SpeciesSummation<SPHBody, DiffusionParticles>> calculate_averaged_(diffusion_body, "Phi");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
	RestartIO	restart_io(io_environment, sph_system.real_bodies_);
	ObservedQuantityRecording<Real> write_solid_("Phi", io_environment, _observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	setup_diffusion_initial_condition.exec();
	setup_boundary_condition.exec();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
		diffusion_body.updateCellLinkedList();
		diffusion_body_complex.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = sph_system.RestartStep();
	Real T0 = 10;
	Real End_Time = T0;
	int restart_output_interval = 1000;
	Real dt = 0.0;
	Real current_averaged_ = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	std::string filefullpath_nonopt_ = io_environment.output_folder_ + "/" + "nonopt_.dat";
	std::ofstream out_file_nonopt_(filefullpath_nonopt_.c_str(), std::ios::app);

	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		dt = get_time_step_size.exec();
		if (ite % 500 == 0)
		{
			write_states.writeToFile(ite);
			write_solid_.writeToFile(ite);

			current_averaged_ = calculate_averaged_.exec();
			out_file_nonopt_ << std::fixed << std::setprecision(12) << ite << "   " << current_averaged_ << "\n";

			std::cout << "N= " << ite << " Time: " << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
			std::cout << "The averaged  is " << calculate_averaged_.exec() << std::endl;
		}

		_splitting.exec(dt);
		ite++; GlobalStaticVariables::physical_time_ += dt;

		if (ite % restart_output_interval == 0)
		{
			restart_io.writeToFile(ite);
		}
	}
	TickCount t4 = TickCount::now();
	TickCount::interval_t tt;
	tt = t4 - t1;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	std::cout << "Total physical time for computation: " << GlobalStaticVariables::physical_time_ << " seconds." << std::endl;

	EXPECT_NEAR(619.124, calculate_averaged_.exec(), 0.01);
}

int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
