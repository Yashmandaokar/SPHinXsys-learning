#include "surface_tension.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
SurfaceTensionStress::
    SurfaceTensionStress(BaseContactRelation &contact_relation, StdVec<Real> contact_surface_tension)
    : LocalDynamics(contact_relation.getSPHBody()), FluidContactData(contact_relation)
{
    particles_->registerVariable(color_gradient_, "ColorGradient");
    particles_->registerSortableVariable<Vecd>("ColorGradient");
    particles_->addVariableToWrite<Vecd>("ColorGradient");
    particles_->registerVariable(surface_tension_stress_, "SurfaceTensionStress");
    particles_->registerSortableVariable<Matd>("SurfaceTensionStress");
    particles_->addVariableToWrite<Matd>("SurfaceTensionStress");
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_surface_tension_.push_back(contact_surface_tension[k]);
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
    }
}
//=================================================================================================//
void SurfaceTensionStress::interaction(size_t index_i, Real dt)
{
    color_gradient_[index_i] = ZeroData<Vecd>::value;
    surface_tension_stress_[index_i] = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd weighted_color_gradient = ZeroData<Vecd>::value;
        Real contact_fraction_k = contact_fraction_[k];
        Real surface_tension_k = contact_surface_tension_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            weighted_color_gradient -= contact_fraction_k *
                                       contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
        }
        color_gradient_[index_i] = weighted_color_gradient;
        Real norm = weighted_color_gradient.norm();
        surface_tension_stress_[index_i] += surface_tension_k / (norm + Eps) *
                                            (norm * norm * Matd::Identity() -
                                             weighted_color_gradient * weighted_color_gradient.transpose());
    }
}
//=================================================================================================//
SurfaceStressAcceleration<Inner<>>::SurfaceStressAcceleration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      rho_(particles_->rho_), mass_(particles_->mass_), force_prior_(particles_->force_prior_),
      color_gradient_(*particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_tension_stress_(*particles_->getVariableByName<Matd>("SurfaceTensionStress")) {}
//=================================================================================================//
void SurfaceStressAcceleration<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        summation += mass_[index_i] * inner_neighborhood.dW_ijV_j_[n] *
                     (surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
                     inner_neighborhood.e_ij_[n];
    }
    force_prior_[index_i] += summation / rho_[index_i];
}
//=================================================================================================//
SurfaceStressAcceleration<Contact<>>::SurfaceStressAcceleration(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), FluidContactData(contact_relation),
      rho_(particles_->rho_), mass_(particles_->mass_), force_prior_(particles_->force_prior_),
      color_gradient_(*particles_->getVariableByName<Vecd>("ColorGradient")),
      surface_tension_stress_(*particles_->getVariableByName<Matd>("SurfaceTensionStress"))
{
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
        contact_color_gradient_.push_back(
            contact_particles_[k]->getVariableByName<Vecd>("ColorGradient"));
        contact_surface_tension_stress_.push_back(
            contact_particles_[k]->getVariableByName<Matd>("SurfaceTensionStress"));
    }
}
//=================================================================================================//
void SurfaceStressAcceleration<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real contact_fraction_k = contact_fraction_[k];
        StdLargeVec<Vecd> &contact_color_gradient_k = *(contact_color_gradient_[k]);
        StdLargeVec<Matd> &contact_surface_tension_stress_k = *(contact_surface_tension_stress_[k]);
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real mismatch = 1.0 - 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]).dot(e_ij) * r_ij;
            summation += mass_[index_i] * contact_neighborhood.dW_ijV_j_[n] *
                         (-0.1 * mismatch * Matd::Identity() +
                          (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] +
                          contact_surface_tension_stress_k[index_j] * contact_fraction_k) *
                         contact_neighborhood.e_ij_[n];
        }
    }
    force_prior_[index_i] += summation / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
