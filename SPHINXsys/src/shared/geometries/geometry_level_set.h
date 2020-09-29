/**
* @file geometry_level_set.h
* @brief Here, we define geometry based on level set technique.
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/

#pragma once

#include "geometry.h"

#include <string>

namespace SPH {

	class SPHBody;
	class BaseLevelSet;
	/**
	 * @class LevelSetComplexShape
	 * @brief the final geomtrical definition of the SPHBody based on a narrow band level set function
	 * generated from the original ComplexShape
	 */
	class LevelSetComplexShape : public ComplexShape
	{
	public:
		LevelSetComplexShape(SPHBody* sph_body, ComplexShape &complex_shape, bool isCleaned = false);
		virtual ~LevelSetComplexShape() {};

		virtual bool checkContain(Vecd input_pnt, bool BOUNDARY_INCLUDED = true) override;
		virtual bool checkNotFar(Vecd input_pnt, Real threshold) override;
		virtual Vecd findClosestPoint(Vecd input_pnt) override;
		virtual Real findSignedDistance(Vecd input_pnt) override;
		virtual Vecd findNormalDirection(Vecd input_pnt) override;
		virtual Vecd computeKernelIntegral(Vecd input_pnt, Kernel * kernel) override;
	protected:
		BaseLevelSet* level_set_;	/**< narrow bounded levelset mesh. */
	};
}

