#ifndef RANS_DYNAMICS_H
#define RANS_DYNAMICS_H
#include "eulerian_fluid_integration.h"
#include "all_particle_dynamics.h"
#include "riemann_solver.h"
#include "fluid_integration.hpp"

#include "base_fluid_dynamics.h"
#include "force_prior.h"
#include "viscous_dynamics.h"
namespace SPH
{
namespace fluid_dynamics
{ 
template <typename... InteractionTypes>
class TurbulentViscousForceInFVM;

    template <class DataDelegationType>
class TurbulentViscousForceInFVM<DataDelegationType>
        : public LocalDynamics, public DataDelegationType
    {
      public:
        template <class BaseRelationType>
        explicit TurbulentViscousForceInFVM(BaseRelationType &base_relation);
        virtual ~TurbulentViscousForceInFVM(){};

      protected:
        StdLargeVec<Real> &rho_, &mass_, &Vol_, &mu_tprof_;
        StdLargeVec<Vecd> &vel_, &turbulent_viscous_force_;
        Real smoothing_length_;
};

    template <>
    class TurbulentViscousForceInFVM<Inner<>>
    : public TurbulentViscousForceInFVM<DataDelegateInner>, public ForcePrior
    {
      public:
        explicit TurbulentViscousForceInFVM(BaseInnerRelation &inner_relation);
        virtual ~TurbulentViscousForceInFVM(){};
        void interaction(size_t index_i, Real dt = 0.0);
    };
    using TurbulentViscousForceInner = TurbulentViscousForceInFVM<Inner<>>;
    //=================================================================================================//
    template <typename... InteractionTypes>
    class TkeGradientForceInFVM;

    template <class DataDelegationType>
    class TkeGradientForceInFVM<DataDelegationType>
        : public LocalDynamics, public DataDelegationType
    {
      public:
        template <class BaseRelationType>
        explicit TkeGradientForceInFVM(BaseRelationType &base_relation);
        virtual ~TkeGradientForceInFVM(){};

      protected:
        StdLargeVec<Real> &rho_, &mass_, &Vol_, &Kprof_;
        StdLargeVec<Vecd> &tke_gradient_force_;
    };

    template <>
    class TkeGradientForceInFVM<Inner<>>
        : public TkeGradientForceInFVM<DataDelegateInner>, public ForcePrior
    {
      public:
        explicit TkeGradientForceInFVM(BaseInnerRelation &inner_relation);
        virtual ~TkeGradientForceInFVM(){};
        void interaction(size_t index_i, Real dt = 0.0);
    };
    using TkeGradientForceInner = TkeGradientForceInFVM<Inner<>>;


    //=================================================================================================//
    
    /*   
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
          StdLargeVec<Real> &Kprof_, &Epsprof_, &mu_tprof_, &walladjacentcellflag_, &Tau_wall_, &wallfacearea_;
            StdLargeVec<Vecd> &meanvel_advection_, &viscous_dissipation_, &pressuregrad_, &tkegrad_, &shear_force_;
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
       */
        //=================================================================================================//
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // RANS_DYNAMICS_H
