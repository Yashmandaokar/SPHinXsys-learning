/**
* @file 	base_mesh.hpp
* @brief 	This is the implementation of the template function and class for base mesh
* @author	Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once

#include "base_mesh.h"
//=================================================================================================//
namespace SPH {
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<class DataType>
	DataType BaseDataPackage<PKG_SIZE, ADDRS_SIZE>
		::probeDataPackage(PackageDataAddress<DataType>& pkg_data_addrs, Vecd& position)
	{
		Vecu grid_idx = GridIndexFromPosition(position);
		Vecd grid_pos = GridPositionFromIndex(grid_idx);
		Vecd alpha = (position - grid_pos) / grid_spacing_;
		Vecd beta = Vec2d(1.0) - alpha;

		DataType bilinear
			= *pkg_data_addrs[grid_idx[0]][grid_idx[1]] * beta[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1]] * alpha[0] * beta[1]
			+ *pkg_data_addrs[grid_idx[0]][grid_idx[1] + 1] * beta[0] * alpha[1]
			+ *pkg_data_addrs[grid_idx[0] + 1][grid_idx[1] + 1] * alpha[0] * alpha[1];

		return  bilinear;
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename InDataType, typename OutDataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeGradient(PackageDataAddress<InDataType>& in_pkg_data_addrs,
			PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
			{
				Real dphidx = (*in_pkg_data_addrs[i + 1][j] - *in_pkg_data_addrs[i - 1][j]);
				Real dphidy = (*in_pkg_data_addrs[i][j + 1] - *in_pkg_data_addrs[i][j - 1]);
				*out_pkg_data_addrs[i][j] = Vecd(dphidx, dphidy);
			}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename InDataType, typename OutDataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		computeNormalizedGradient(PackageDataAddress<InDataType>& in_pkg_data_addrs,
			PackageDataAddress<OutDataType> out_pkg_data_addrs, Real dt)
	{
		for (int i = 1; i != PKG_SIZE + 1; ++i)
			for (int j = 1; j != PKG_SIZE + 1; ++j)
			{
				Real dphidx = (*in_pkg_data_addrs[i + 1][j] - *in_pkg_data_addrs[i - 1][j]);
				Real dphidy = (*in_pkg_data_addrs[i][j + 1] - *in_pkg_data_addrs[i][j - 1]);
				Real norm = sqrt(dphidx * dphidx + dphidy * dphidy) + TinyReal;
				*out_pkg_data_addrs[i][j] = Vecd(dphidx, dphidy) / norm;
			}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		initializePackageDataAddress(PackageData<DataType>& pkg_data,
			PackageDataAddress<DataType>& pkg_data_addrs)
	{
		for (int i = 0; i != ADDRS_SIZE; ++i)
			for (int j = 0; j != ADDRS_SIZE; ++j)
			{
				pkg_data_addrs[i][j] = &pkg_data[0][0];
			}
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	DataType  BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		CornerAverage(PackageDataAddress<DataType>& pkg_data_addrs, Veci addrs_index, Veci corner_direction)
	{
		DataType average(0);
		for (int i = 0; i != 2; ++i)
			for (int j = 0; j != 2; ++j)
			{
				int x_index = addrs_index[0] + i * corner_direction[0];
				int y_index = addrs_index[1] + j * corner_direction[1];
				average += *pkg_data_addrs[x_index][y_index];
			}
		return average * 0.25;
	}
	//=================================================================================================//
	template<int PKG_SIZE, int ADDRS_SIZE>
	template<typename DataType>
	void BaseDataPackage<PKG_SIZE, ADDRS_SIZE>::
		assignPackageDataAddress(PackageDataAddress<DataType>& pkg_data_addrs, Vecu& addrs_index,
			PackageData<DataType>& pkg_data, Vecu& data_index)
	{
		pkg_data_addrs[addrs_index[0]][addrs_index[1]] = &pkg_data[data_index[0]][data_index[1]];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> 
	template<typename DataType, typename PackageDataType, PackageDataType DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::
		DataValueFromGlobalIndex(Vecu global_data_index)
	{
		Vecu pkg_index_(0);
		Vecu local_data_index(0);
		for (int n = 0; n != 2; n++)
		{
			size_t cell_index_in_this_direction = global_data_index[n] / pkg_size_;
			pkg_index_[n] = cell_index_in_this_direction;
			local_data_index[n] = global_data_index[n] - cell_index_in_this_direction * pkg_size_;
		}
		PackageDataType& data = data_pkg_addrs_[pkg_index_[0]][pkg_index_[1]]->*MemPtr;
		return data[local_data_index[0]][local_data_index[1]];
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType> 
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::initializePackageAddressesInACell(Vecu cell_index)
	{
		int i = (int)cell_index[0];
		int j = (int)cell_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		if (data_pkg->is_inner_pkg_) {
			for (int l = 0; l != pkg_addrs_size_; ++l)
				for (int m = 0; m != pkg_addrs_size_; ++m) {
					pair<int, int>  x_pair = CellShiftAndDataIndex(l);
					pair<int, int>  y_pair = CellShiftAndDataIndex(m);
					data_pkg->assignAllPackageDataAddress(Vecu(l, m), 
						data_pkg_addrs_[i + x_pair.first][j + y_pair.first],
						Vecu(x_pair.second, y_pair.second));
				}
		}
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::allocateMeshDataMatrix()
	{
		Allocate2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	void MeshWithDataPackages<BaseMeshType, DataPackageType>::deleteMeshDataMatrix()
	{
		Delete2dArray(data_pkg_addrs_, BaseMeshType::number_of_cells_);
	}
	//=================================================================================================//
	template<class BaseMeshType, class DataPackageType>
	template<class DataType, typename PackageDataAddressType, PackageDataAddressType DataPackageType:: * MemPtr>
	DataType MeshWithDataPackages<BaseMeshType, DataPackageType>::probeMesh(Vecd& position)
	{
		Vecu grid_index = BaseMeshType::GridIndexFromPosition(position);
		size_t i = grid_index[0];
		size_t j = grid_index[1];

		DataPackageType* data_pkg = data_pkg_addrs_[i][j];
		PackageDataAddressType& pkg_data_addrs = data_pkg->*MemPtr;
		return data_pkg->is_inner_pkg_ ?
			data_pkg->DataPackageType::template probeDataPackage<DataType>(pkg_data_addrs, position)
			: *pkg_data_addrs[0][0];
	}
	//=================================================================================================//
}
//=================================================================================================//
