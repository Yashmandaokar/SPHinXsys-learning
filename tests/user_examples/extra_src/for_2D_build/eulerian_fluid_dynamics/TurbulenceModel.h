#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H
#include "eulerian_fluid_dynamics.h"
#include "common_shared_FVM_classes.h"

namespace SPH
{
    namespace fluid_dynamics
    {

        class MeanVelocity : public BaseIntegration
        {
        public:
            explicit MeanVelocity(BaseInnerRelation& inner_relation);
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd> meanvelocity_;
            StdLargeVec<Real>& Vol_;

        };

        class BaseTurbulence :  public MeanVelocity 
        {
            public:
            explicit BaseTurbulence(BaseInnerRelation &inner_relation);
            virtual ~BaseTurbulence() {};

            protected:
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_;
            StdLargeVec<Vecd> &meanvelocity_;
            StdLargeVec<Real> &mu_t_, &Vol_, &K_, &Eps_;
        };

        class WeaklyCompressibleFluidInitialCondition_KEps : public FluidInitialCondition
        {
            public:
            explicit WeaklyCompressibleFluidInitialCondition_KEps(SPHBody &sph_body)
            : FluidInitialCondition(sph_body), rho_(particles_->rho_),
            p_(*particles_->getVariableByName<Real>("Pressure"))
            {};

            protected:
            StdLargeVec<Real> &rho_, &p_;
        };

        
        class KEpsilonStd1stHalf : public BaseTurbulence
        {
            public:
            explicit KEpsilonStd1stHalf(BaseInnerRelation &inner_relation);
            virtual ~KEpsilonStd1stHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            StdLargeVec<Real> dK_dt_;
        };

        class KEpsilonStd2ndHalf : public BaseTurbulence
        {
            public:
            explicit KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation);
            virtual ~KEpsilonStd2ndHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            StdLargeVec<Real> dEps_dt_;
        };

        class NoRiemannSolverInCompressobleRANS
        {
            public:
            explicit NoRiemannSolverInCompressobleRANS(Fluid& fluid_i, Fluid& fluid_j, Real limiter_parameter = 15.0)
                : fluid_i_(fluid_i), fluid_j_(fluid_j), limiter_parameter_(limiter_parameter) {};
            FluidStarState getInterfaceState(const FluidState& state_i, const FluidState& state_j, const Vecd& e_ij);
            
            protected:
            Fluid& fluid_i_, & fluid_j_;
            Real limiter_parameter_;
        };
        

        template <class RiemannSolverType>
        class EulerianIntegration1stHalfRANS : public BaseTurbulence
        {
            public:
            explicit EulerianIntegration1stHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
            virtual ~EulerianIntegration1stHalfRANS(){};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real limiter_input_;
            RiemannSolverType riemann_solver_;
            StdLargeVec<Vecd> &acc_prior_;
            StdLargeVec<Vecd> mom_, dmom_dt_;
        };
        /* define the mostly used pressure relaxation scheme using Riemann solver */
        using EulerianIntegration1stHalfRANSNoRiemannSolver = EulerianIntegration1stHalfRANS<NoRiemannSolverInCompressobleRANS>;

        /*
         * @class EulerianIntegration2ndHalf
         * @brief  Template density relaxation scheme with different Riemann solver*/
         
        template <class RiemannSolverType>
        class EulerianIntegration2ndHalfRANS : public BaseTurbulence
        {
        public:
            explicit EulerianIntegration2ndHalfRANS(BaseInnerRelation &inner_relation, Real limiter_parameter = 15.0);
            virtual ~EulerianIntegration2ndHalfRANS(){};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            Real limiter_input_;
            RiemannSolverType riemann_solver_;
        };
        using EulerianIntegration2ndHalfRANSNoRiemannSolver = EulerianIntegration2ndHalfRANS<NoRiemannSolverInCompressobleRANS>;

       

       class StdWallFunction : public BaseTurbulence
       {
            public:
            explicit StdWallFunction (BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
            virtual ~StdWallFunction(){};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);
            

            protected:
            GhostCreationFromMesh &ghost_creator_;
            Real vonkar_, E_;
            StdLargeVec<Real> dK_dt_, Esp_p_;
        };

    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCEMODEL_H
