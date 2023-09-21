#include "fluid_body.h"

#include "base_material.h"
#include "base_particles.hpp"
#include "cell_linked_list.h"

namespace SPH
{
//=================================================================================================//
FluidBody::FluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
    : RealBody(system, shape_ptr) {}
//=================================================================================================//

/**
	 * @class EulerianFluidBody
	 * @brief Eulerian Fluid body uses smoothing length to particle spacing 1.3 
	 */
	EulerianFluidBody::EulerianFluidBody(SPHSystem &system, SharedPtr<Shape> shape_ptr)
		: FluidBody(system, shape_ptr)
	{
		defineAdaptation<SPHAdaptation>(1.3);
	}
} // namespace SPH
