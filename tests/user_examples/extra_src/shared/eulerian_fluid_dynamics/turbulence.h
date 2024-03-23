
#ifndef TURBULENCE_H
#define TURBULENCE_H


#include "eulerian_fluid_dynamics.h"

namespace SPH
{
    namespace fluid_dynamics
    {
        /*
        @class WeaklyCompressibleFluidInitialCondition
        @brief  Set initial condition for a Eulerian weakly compressible fluid body.
        This is a abstract class to be override for case specific initial conditions}*/
        /*
        class WeaklyCompressibleFluidInitialCondition_KEps : public FluidInitialCondition
        {
            public:
            explicit WeaklyCompressibleFluidInitialCondition_KEps(SPHBody &sph_body)
            : FluidInitialCondition(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), 
            mom_(*particles_->getVariableByName<Vecd>("Momentum")), K_(*particles_->getVariableByName<Vecd>("TKE")), 
            Epsilon_(*particles_->getVariableByName<Vecd>("DissipationRate"))
            {};

            protected:
            StdLargeVec<Real> &rho_, &p_;
            StdLargeVec<Vecd> &mom_, &K_, &Epsilon_;
        };*/
        //=================================================================================================//
        class TurbulenceModel_K_Epsilon 
        {
            public: 
            explicit TurbulenceModel_K_Epsilon(Real& C_mu);
            virtual ~TurbulenceModel_K_Epsilon(){};
            Real TurbulentViscosity(Real& rho_star, Vecd& K_initial, Vecd& Epsilon_initial);
            Matd Production_TKE(Vecd& K_initial, FluidStarState& interface_state, Real& dW_ijV_j, Real& mu_t_initial);
            Vecd TkeTransportEqn(Vecd& K_initial, Vecd& Epsilon_initial, Real& rho_star, FluidStarState& interface_state,
                                Real& mu_t_initial, Real& mu_, Vecd& TKE_change_rate, Real& dW_ijV_j, Real &r_ij, Vecd &e_ij, Real& dt, Matd& P_TKE);
            Vecd EpsilonTransportEqn(Vecd& K_updated, Vecd& Epsilon_initial, Real& rho_star, FluidStarState& interface_state,
                    Real& mu_t_initial, Real& mu_, Vecd& Epsilon_change_rate, Real& dW_ijV_j, Real &r_ij, Vecd &e_ij, Real& dt, Matd& P_TKE);
            Real CorrectedTurbulentViscosity(Real& rho_star, Vecd& K_updated, Vecd& Epsilon_updated);
            Matd ReynoldStressTensor(Vecd& K_updated, FluidStarState& interface_state, Real& dW_ijV_j, Real& mu_t_corrected);

            protected:
            Real &C_mu_;
           
        };
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCE_H
