#ifndef RANS_FLUID_INTEGRATION_H
#define RANS_FLUID_INTEGRATION_H
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
        class EulerianIntegrationRANS : public BaseIntegration<DataDelegationType>
        {
        public:
            template <class BaseRelationType>
            explicit EulerianIntegrationRANS(BaseRelationType& base_relation);
            virtual ~EulerianIntegrationRANS() {};

        protected:
            StdLargeVec<Vecd> &mom_, &dmom_dt_;
            StdLargeVec<Real> &dmass_dt_, &Vol_;
            Real Cmu_;
            //StdLargeVec<Real> &K_, &Eps_, &mu_t_;
        };


        template <typename... InteractionTypes>
        class EulerianIntegration1stHalfRANS;

        template <class RiemannSolverType>
        class EulerianIntegration1stHalfRANS<Inner<>, RiemannSolverType>
            : public EulerianIntegrationRANS<DataDelegateInner>
        {
        public:
            explicit EulerianIntegration1stHalfRANS(BaseInnerRelation& inner_relation, Real limiter_parameter = 15.0);
            template <typename BodyRelationType, typename FirstArg>
            explicit EulerianIntegration1stHalfRANS(ConstructorArgs<BodyRelationType, FirstArg> parameters)
                : EulerianIntegration1stHalfRANS(parameters.body_relation_, std::get<0>(parameters.others_)) {};
            virtual ~EulerianIntegration1stHalfRANS() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            RiemannSolverType riemann_solver_;
            StdLargeVec<Real> &Kprof_, &Epsprof_, &mu_tprof_;
            StdLargeVec<Vecd> &vel_prof_, &meanvel_advection_, &viscous_dissipation_, &pressuregrad_, &tkegrad_;
            StdLargeVec<Matd> &vel_gradient_mat_;
        };
        using EulerianIntegration1stHalfInnerRiemannRANS = EulerianIntegration1stHalfRANS<Inner<>, AcousticRiemannSolver>;

        using BaseEulerianIntegrationWithWallRANS = InteractionWithWall<EulerianIntegrationRANS>;
        template <class RiemannSolverType>
        class EulerianIntegration1stHalfRANS<Contact<Wall>, RiemannSolverType>
            : public BaseEulerianIntegrationWithWallRANS
        {
        public:
            EulerianIntegration1stHalfRANS(BaseContactRelation& wall_contact_relation, Real limiter_parameter = 15.0);
            template <typename BodyRelationType, typename FirstArg>
            explicit EulerianIntegration1stHalfRANS(ConstructorArgs<BodyRelationType, FirstArg> parameters)
                : EulerianIntegration1stHalfRANS(parameters.body_relation_, std::get<0>(parameters.others_)) {};
            virtual ~EulerianIntegration1stHalfRANS() {};
            void interaction(size_t index_i, Real dt = 0.0);

        protected:
            RiemannSolverType riemann_solver_;
            StdLargeVec<Real> &K_, &Eps_, &mu_t_;
        };
        using EulerianIntegration1stHalfWithWallRiemannRANS =
            ComplexInteraction<EulerianIntegration1stHalfRANS<Inner<>, Contact<Wall>>, AcousticRiemannSolver>;

        //=================================================================================================//
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // RANS_FLUID_INTEGRATION_H
