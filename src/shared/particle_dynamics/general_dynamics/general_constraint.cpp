#include "general_constraint.h"

namespace SPH
{
//=================================================================================================//
ShapeSurfaceBounding::ShapeSurfaceBounding(NearShapeSurface &near_shape_surface)
    : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
      GeneralDataDelegateSimple(sph_body_), pos_(particles_->pos_),
      constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
{
    level_set_shape_ = &near_shape_surface.level_set_shape_;
}
//=================================================================================================//
void ShapeSurfaceBounding::update(size_t index_i, Real dt)
{
    Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

    if (phi > -constrained_distance_)
    {
        Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
        pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
    }
}
//=================================================================================================//
} // namespace SPH
