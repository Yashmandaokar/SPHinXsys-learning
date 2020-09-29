/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file base_particle_dynamics.h
* @brief This is for the base classes of particle dynamics, which describe the
* interaction between particles. These interactions are used to define  
* differential operators for surface forces or fluxes in continuum mechanics
* @author  Xiangyu Hu, Luhui Han and Chi Zhang
*/
#pragma once
#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_particles.h"
#include "all_materials.h"
#include "neighbor_relation.h"
#include "all_types_of_bodies.h"
#include "all_meshes.h"
#include "external_force.h"
#include "body_relation.h"
#include <functional>

using namespace std::placeholders;

namespace SPH 
{
	/** Functor for operation of inner particles. */
	typedef std::function<void(size_t, Real)> InnerFunctor;
	/** Functors for reducing operation of inner particles. */
	template <class ReturnType>
	using ReduceFunctor = std::function<ReturnType(size_t, Real)>;

	/** Iterators for inner functors. sequential computing. */
	void InnerIterator(size_t number_of_particles, InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors. parallel computing. */
	void InnerIterator_parallel(size_t number_of_particles, InnerFunctor &inner_functor, Real dt = 0.0);

	/** Iterators for reduce functors. sequential computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt = 0.0);
	/** Iterators for reduce functors. parallel computing. */
	template <class ReturnType, typename ReduceOperation>
	ReturnType ReduceIterator_parallel(size_t number_of_particles, ReturnType temp,
		ReduceFunctor<ReturnType> &reduce_functor, ReduceOperation &reduce_operation, Real dt = 0.0);

	/** Functor for configuration operation. */
	typedef std::function<void(CellList*, Real)> CellListFunctor;
	/** Iterators for inner functors with splitting for configuration dynamics. sequential computing. */
	void CellListIteratorSplitting(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting for configuration dynamics. parallel computing. */
	void CellListIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		CellListFunctor& cell_list_functor, Real dt = 0.0);

	/** Iterators for inner functors with splitting. sequential computing. */
	void InnerIteratorSplitting(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. parallel computing. */
	void InnerIteratorSplitting_parallel(SplitCellLists& split_cell_lists,
		InnerFunctor &inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. sequential computing. */
	void InnerIteratorSplittingSweeping(SplitCellLists& split_cell_lists,
		InnerFunctor& inner_functor, Real dt = 0.0);
	/** Iterators for inner functors with splitting. parallel computing. */
	void InnerIteratorSplittingSweeping_parallel(SplitCellLists& split_cell_lists,
		InnerFunctor& inner_functor, Real dt = 0.0);


	/** A Functor for Summation */
	template <class ReturnType>
	struct ReduceSum { ReturnType operator () (ReturnType x, ReturnType y) const { return x + y; }; };
	/** A Functor for Maximum */
	struct ReduceMax { Real operator () (Real x, Real y) const { return SMAX(x, y); }; };
	/** A Functor for Minimum */
	struct ReduceMin { Real operator () (Real x, Real y) const { return SMIN(x, y); }; };
	/** A Functor for OR operator */
	struct ReduceOR { bool operator () (bool x, bool y) const { return x || y; }; };
	/** A Functor for lower bound */
	struct ReduceLowerBound {
		Vecd operator () (Vecd x, Vecd y) const {
			Vecd lower_bound;
			for (int i = 0; i < lower_bound.size(); ++i) lower_bound[i] = SMIN(x[i], y[i]);
			return lower_bound;
		}; 
	};
	/** A Functor for upper bound */
	struct ReduceUpperBound {
		Vecd operator () (Vecd x, Vecd y) const {
			Vecd upper_bound;
			for (int i = 0; i < upper_bound.size(); ++i) upper_bound[i] = SMAX(x[i], y[i]);
			return upper_bound;
		};
	};

	/**
	 * @class GlobalStaticVariables
	 * @brief A place to put all global variables
	 */
	class GlobalStaticVariables
	{
	public:
		explicit GlobalStaticVariables() {};
		virtual ~GlobalStaticVariables() {};

		/** the physical time is global value for all dynamics */
		static Real physical_time_;
	};

	/**
	* @class ParticleDynamics
	* @brief The base class for all particle dynamics
	* This class contains the only two interface functions available
	* for particle dynamics. An specific implementation should be realized.
	*/
	template <class ReturnType = void>
	class ParticleDynamics : public GlobalStaticVariables
	{
	public:
		/** Constructor */
		explicit ParticleDynamics(SPHBody* sph_body) 
			: GlobalStaticVariables(), sph_body_(sph_body),
			split_cell_lists_(sph_body->split_cell_lists_),
			mesh_cell_linked_list_(sph_body->mesh_cell_linked_list_) {};
		virtual ~ParticleDynamics() {};

		/** The only two functions can be called from outside
		  * One is for sequential execution, the other is for parallel. */
		virtual ReturnType exec(Real dt = 0.0) = 0;
		virtual ReturnType parallel_exec(Real dt = 0.0) = 0;
	protected:
		SPHBody* sph_body_;
		SplitCellLists& split_cell_lists_;
		BaseMeshCellLinkedList* mesh_cell_linked_list_;

		void setBodyUpdated() { sph_body_->setNewlyUpdated(); };
		/** the function for set global parameters for the particle dynamics */
		virtual void setupDynamics(Real dt = 0.0) {};
	};

	/**
	* @class DataDelegateSimple
	* @brief prepare data for simple particle dynamics.
	*/
	template <class BodyType = SPHBody, 
			  class ParticlesType = BaseParticles, 
			  class MaterialType = BaseMaterial>
	class DataDelegateSimple
	{
	public:
		/** Constructor */
		explicit DataDelegateSimple(SPHBody* body) :
			body_(dynamic_cast<BodyType*>(body)),
			particles_(dynamic_cast<ParticlesType*>(body->base_particles_)),
			material_(dynamic_cast<MaterialType*>(body->base_particles_->base_material_)) {};
		virtual ~DataDelegateSimple() {};

	protected:
		BodyType* body_;
		ParticlesType* particles_;
		MaterialType* material_;
	};

	/**
	* @class DataDelegateInner
	* @brief prepare data for inner particle dynamics
	*/
	template <class BodyType = SPHBody,
			  class ParticlesType = BaseParticles,
			  class MaterialType = BaseMaterial>
	class DataDelegateInner
	{
	public:
		explicit DataDelegateInner(SPHBodyInnerRelation* body_inner_relation) : 
			body_(dynamic_cast<BodyType*>(body_inner_relation->sph_body_)),
			particles_(dynamic_cast<ParticlesType*>(body_->base_particles_)),
			material_(dynamic_cast<MaterialType*>(body_->base_particles_->base_material_)),
			inner_configuration_(body_inner_relation->inner_configuration_) {};
		virtual ~DataDelegateInner() {};
	protected:
		BodyType* body_;
		ParticlesType* particles_;
		MaterialType* material_;
		/** inner configuration of the designated body */
		ParticleConfiguration& inner_configuration_;
	};

	/**
	* @class DataDelegateContact
	* @brief prepare data for contact particle dynamics
	*/
	template <class BodyType = SPHBody,
			  class ParticlesType = BaseParticles,
			  class MaterialType = BaseMaterial,
			  class ContactBodyType = SPHBody,
			  class ContactParticlesType = BaseParticles,
			  class ContactMaterialType = BaseMaterial>
	class DataDelegateContact
	{
	public:
		explicit DataDelegateContact(SPHBodyContactRelation* body_contact_relation);
		virtual ~DataDelegateContact() {};
	protected:
		BodyType* body_;
		ParticlesType* particles_;
		MaterialType* material_;

		StdVec<ContactBodyType*>  contact_bodies_;
		StdVec<ContactParticlesType*>  contact_particles_;
		StdVec<ContactMaterialType*>  contact_material_;
		/** Configurations for particle interaction between bodies. */
		ContatcParticleConfiguration& contact_configuration_;
	};

	/**
	* @class DataDelegateComplex
	* @brief prepare data for complex particle dynamics
	*/
	template <class BodyType = SPHBody,
		class ParticlesType = BaseParticles,
		class MaterialType = BaseMaterial,
		class ContactBodyType = SPHBody,
		class ContactParticlesType = BaseParticles,
		class ContactMaterialType = BaseMaterial>
	class DataDelegateComplex
	{
	public:
		explicit DataDelegateComplex(SPHBodyComplexRelation* body_complex_relation);
		virtual ~DataDelegateComplex() {};
	protected:
		BodyType* body_;
		ParticlesType* particles_;
		MaterialType* material_;
		/** inner configuration of the designated body */
		ParticleConfiguration& inner_configuration_;

		StdVec<ContactBodyType*>  contact_bodies_;
		StdVec<ContactParticlesType*>  contact_particles_;
		StdVec<ContactMaterialType*>  contact_material_;
		/** Configurations for particle interaction between bodies. */
		ContatcParticleConfiguration& contact_configuration_;
	};
}
