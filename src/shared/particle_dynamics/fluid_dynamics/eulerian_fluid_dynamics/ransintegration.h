#ifndef RANSINTEGRATION_H
#define RANSINTEGRATION_H
#include "eulerian_fluid_integration.h"
#include "all_particle_dynamics.h"
#include "riemann_solver.h"
#include "fluid_integration.hpp"
namespace SPH
{
    namespace fluid_dynamics
    {   
        //=================================================================================================//
        template <class DataDelegationType>
        class IntegrationRANS : public BaseIntegration<DataDelegationType>
        {
        public:
            template <class BaseRelationType>
            explicit IntegrationRANS(BaseRelationType& base_relation);
            virtual ~IntegrationRANS() {};

        protected:
            StdLargeVec<Vecd> &mom_, &dmom_dt_;
            StdLargeVec<Real> &dmass_dt_, &Vol_;
        };


        template <typename... InteractionTypes>
        class Integration1stHalfRANS;

        template <class RiemannSolverType>
        class Integration1stHalfRANS<Inner<>, RiemannSolverType>
            : public IntegrationRANS<DataDelegateInner>
        {
        public:
            explicit Integration1stHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
            template <typename BodyRelationType, typename FirstArg>
            explicit Integration1stHalfRANS(ConstructorArgs<BodyRelationType, FirstArg> parameters)
                : Integration1stHalfRANS(parameters.body_relation_, std::get<0>(parameters.others_)){};
            virtual ~Integration1stHalfRANS(){};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            RiemannSolverType riemann_solver_;
        };
        using Integration1stHalfInnerRiemannRANS = Integration1stHalfRANS<Inner<>, AcousticRiemannSolver>;
        //=================================================================================================//
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // RANSINTEGRATION_H
