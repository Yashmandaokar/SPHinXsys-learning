#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H
#include "fluid_integration.hpp"
#include "unstructured_mesh.h"
#include "all_particle_dynamics.h"
#include "riemann_solver.h"
#include "rans_fluid_integration.h"

namespace SPH
{
    namespace fluid_dynamics
    {
        class BaseTurbulence : public BaseIntegration<FluidDataInner>
        {
        public:
            explicit BaseTurbulence(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~BaseTurbulence() {};
            

        protected:
            StdLargeVec<Vecd>& mom_, & dmom_dt_;
            StdLargeVec<Real>& dmass_dt_;
            StdLargeVec<Real> &K_prod_p_, &K_prodsum_, &Eps_p_, &Eps_sum_;
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_;
            StdLargeVec<Real> &K_, &Eps_, &mu_t_;
            GhostCreationFromMesh& ghost_creator_;
            
        };
        //=================================================================================================//
        class WallAdjacentCells : public BaseTurbulence
        {
        public:
            explicit WallAdjacentCells(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~WallAdjacentCells() {};
            

        protected:
            StdLargeVec<Real> walladjacentindex_, walladjacentcellflag_, wallghostindex_;
            StdLargeVec<Vecd> walleij_;

            void walladjacentcells();
        };
        //=================================================================================================//
        /*template <>
/*
        template <typename... InteractionTypes>
        class BaseTurbulence;

        template <class DataDelegationType>
        class BaseTurbulence<DataDelegationType>
            : public BaseIntegration<DataDelegationType>
        {
        public:
            template <class BaseRelationType>
            explicit BaseTurbulence(BaseRelationType& base_relation);
            virtual ~BaseTurbulence() {};

        protected:
            StdLargeVec<Vecd> &mom_, &dmom_dt_;
            StdLargeVec<Real> &dmass_dt_, &Vol_;
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_;
            StdLargeVec<Real> &K_, &Eps_, &mu_t_;
        };*/
        //=================================================================================================//
        /*template <>
        class BaseTurbulence<Inner<>>
            : public BaseTurbulence<FluidDataInner>
        {
        public:
            explicit BaseTurbulence(BaseInnerRelation& inner_relation);
            virtual ~BaseTurbulence() {};
        };*/
        //=================================================================================================//
        /*using BaseTurbulenceWithWall = InteractionWithWall<BaseTurbulence>;
        template <>
        class BaseTurbulence<Contact<Wall>> : public BaseTurbulenceWithWall
        {
        public:
            explicit BaseTurbulence(BaseContactRelation& wall_contact_relation);
            virtual ~BaseTurbulence() {};
        };*/
        //=================================================================================================//
        class KEpsilonStd1stHalf : public WallAdjacentCells
        {
            public:
            explicit KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~KEpsilonStd1stHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            StdLargeVec<Real> &dK_dt_;
        };
        //=================================================================================================//
        class KEpsilonStd2ndHalf : public WallAdjacentCells
        {
            public:
            explicit KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~KEpsilonStd2ndHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            StdLargeVec<Real> dEps_dt_;
        };
        //=================================================================================================//
            class StdWallFunctionFVM : public WallAdjacentCells
        {
            public:
            explicit StdWallFunctionFVM(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~StdWallFunctionFVM() {};
            void interaction(size_t index_i, Real dt = 0.0);
           // void update(size_t index_i, Real dt = 0.0);


            protected:
            Real vonkar_, E_;
            //StdLargeVec<Real> &dK_dt_;
        };
        //=================================================================================================//
        using BaseIntegrationFSI = InteractionWithWall<BaseIntegration>;
        class StdWallFunction : public BaseIntegrationFSI//InteractionWithWall<BaseIntegration>
        {
            public:
            StdWallFunction(BaseContactRelation& wall_contact_relation);
            
            virtual ~StdWallFunction() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real vonkar_, E_;
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_;
            StdLargeVec<Real> &K_, &Eps_, &mu_t_, &Esp_p_;
            StdLargeVec<Real> &dK_dt_;
            StdLargeVec<Vecd> &distance_from_wall_;
            SPHBody &bounds_;
            StdVec<StdLargeVec<Vecd>*> wall_pos_;
        };
        //=================================================================================================//
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCEMODEL_H
