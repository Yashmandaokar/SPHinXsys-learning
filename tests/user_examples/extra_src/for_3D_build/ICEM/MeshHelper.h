#ifndef MESHTYPE_H
#define MESHTYPE_H

#include "unstructured_mesh.h"
using namespace std;
namespace SPH
{

    class MeshFileHelpers 
    {
        public:

        MeshFileHelpers()
        {
     
        };
        virtual ~MeshFileHelpers(){};

        static void meshDimension(ifstream& mesh_file, size_t& dimension, string& text_line);
        static void numberofNodes(ifstream& mesh_file, size_t& number_of_points, string& text_line);
        static void nodeCoordinates(ifstream& mesh_file, size_t& node_index, vector<vector<Real>> &point_coordinates_3D_, string& text_line, size_t& dimension, vector<vector<Real>>& point_coordinates);
        static void numberofElements(ifstream& mesh_file, size_t& number_of_elements, string& text_line);
        static void dataStruct(vector<vector<vector<size_t>>>& cell_lists_, vector<vector<size_t>>& elements_nodes_connection_, size_t number_of_elements, size_t mesh_type, size_t dimension, StdLargeVec<Vec3d>& elements_neighbors_connection_);
        static size_t findBoundaryType(string& text_line, size_t boundary_type);
        static Vecd nodeIndex(string& text_line);
        static Vec2d cellIndex(string& text_line);
        static void updateElementsNodesConnection(vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2);
        static void updateCellLists(vector<vector<vector<size_t>>>& cell_lists_, vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type);
        static void updateBoundaryCellLists(vector<vector<vector<size_t>>>& cell_lists_, vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type);
        static void cellCenterCoordinates(vector<vector<size_t>>& elements_nodes_connection_, std::vector<std::vector<long unsigned int>>::size_type& element, vector<vector<Real>>& point_coordinates_3D_, StdLargeVec<Vecd>& elements_center_coordinates_, Vecd& center_coordinate); 
        static void elementVolume(vector<vector<size_t>>& elements_nodes_connection_, std::vector<std::vector<long unsigned int>>::size_type& element, vector<vector<Real>>& point_coordinates_3D_, StdLargeVec<Real>& elements_volumes_);
        static void faceArea(vector<Real>& all_data_of_distance_between_nodes, StdLargeVec<Real>& elements_volumes_, vector<vector<vector<size_t>>>& cell_lists_, vector<vector<Real>>& point_coordinates_3D_);
        static void vtuFileHeader(std::ofstream& out_file);
        static void vtuFileNodeCoordinates(std::ofstream& out_file, vector<vector<Real>>& nodes_coordinates_, vector<vector<size_t>>& elements_nodes_connection_);
        static void vtuFileInformationKey(std::ofstream& out_file);
        static void vtuFileCellConnectivity(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_);
        static void vtuFileOffsets(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_);
        static void vtuFileTypeOfCell(std::ofstream& out_file, vector<vector<size_t>>& elements_nodes_connection_);
    };


}
#endif // MESHTYPE_H