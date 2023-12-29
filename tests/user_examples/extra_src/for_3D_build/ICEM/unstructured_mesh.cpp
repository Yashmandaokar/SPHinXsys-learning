
#include "unstructured_mesh.h"
#include "MeshHelper.h"
#include "MeshFileHelperFluent.h"
namespace SPH
{
    //=================================================================================================//
    void readMeshFile_3d::getDataFromMeshFile3d()
    {
        Real ICEM = 0;
        Real FLUENT = 0;
        ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
        mesh_file.open(full_path_);
        if (mesh_file.fail())
        {
            cout << "Error:Check if the file exists." << endl;
        }
        string text_line;
        /*Check mesh file generation software*/
        (getline(mesh_file, text_line));
        text_line.erase(0, 4);
        istringstream value(text_line);
        if (text_line.find("Created by", 0) != string::npos)
        {
            ICEM = 1;
        }
       
        if (ICEM == 1)
        {

            /*--- Read the dimension of the mesh ---*/
            size_t dimension(0);
            MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

            /*--- Read the node data (index is starting from zero) ---*/
            size_t number_of_points(0);
            MeshFileHelpers::numberofNodes(mesh_file, number_of_points, text_line);

            /* Prepare our data structure for the point coordinates. */
            point_coordinates_3D_.resize(number_of_points);
            for (std::vector<std::vector<double>>::size_type point = 0; point != point_coordinates_3D_.size(); ++point)
            {
                point_coordinates_3D_[point].resize(dimension);
            }

            point_coordinates.resize(dimension);
            for (std::size_t k = 0; k < dimension; k++)
            {
                point_coordinates[k].reserve(number_of_points);
            }

            size_t node_index = 0;
            MeshFileHelpers::nodeCoordinates(mesh_file, node_index, point_coordinates_3D_, text_line, dimension, point_coordinates);

            size_t boundary_type(0);
            size_t number_of_elements(0);
            size_t mesh_type = 4;
            MeshFileHelpers::numberofElements(mesh_file, number_of_elements, text_line);

            /*Preparing and initializing the data structure of cell_lists and element node connection*/
            MeshFileHelpers::dataStruct(cell_lists_, elements_nodes_connection_, number_of_elements, mesh_type, dimension, elements_neighbors_connection_);

            /*--- find the elements lines ---*/
            while (getline(mesh_file, text_line))
            {
                if (text_line.find("(13", 0) != string::npos && text_line.find(")(", 0) != string::npos)
                {
                    boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                    types_of_boundary_condition_.push_back(boundary_type);
                    while (getline(mesh_file, text_line))
                    {
                        if (text_line.find(")", 0) == string::npos)
                        {
                            Vecd nodes = MeshFileHelpers::nodeIndex(text_line);
                            Vec2d cells = MeshFileHelpers::cellIndex(text_line);

                            /*--- build up all topology---*/
                            bool check_neighbor_cell1 = 1;
                            bool check_neighbor_cell2 = 0;
                            for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                            {
                                if (mesh_type == 4)
                                {
                                    if (cells[check_neighbor_cell2] != 0)
                                    {
                                        MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);

                                        /*--- build up all connection data with element and neighbor and nodes---*/
                                        MeshFileHelpers::updateCellLists(cell_lists_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);

                                        MeshFileHelpers::updateBoundaryCellLists(cell_lists_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                    }
                                    if (cells[check_neighbor_cell2] == 0)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                        else
                            break;
                    }
                }
                if (text_line.find("Zone Sections", 0) != string::npos)
                    break;
                if (text_line.find(")") != string::npos)
                    continue;
            }
            cell_lists_.erase(cell_lists_.begin());
        }
        else /*This section is for mesh files created from fluent*/
        {
            size_t dimension(0);
            MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

            /*--- Read the node data (index is starting from zero) ---*/
            size_t number_of_points(0);
            MeshFileHelpersFluent::numberofNodesFluent(mesh_file, number_of_points, text_line);
            point_coordinates_3D_.resize(number_of_points);
            for (std::vector<std::vector<double>>::size_type point = 0; point != point_coordinates_3D_.size(); ++point)
            {
                point_coordinates_3D_[point].resize(dimension);
            }
            point_coordinates.resize(dimension);
            for (std::size_t k = 0; k < dimension; k++)
            {
                point_coordinates[k].reserve(number_of_points);
            }
            
            size_t node_index = 0;
            size_t boundary_type(0);
            size_t number_of_elements(0);
            size_t mesh_type = 4;

            MeshFileHelpersFluent::numberofElementsFluent(mesh_file, number_of_elements, text_line);
            MeshFileHelpers::dataStruct(cell_lists_, elements_nodes_connection_, number_of_elements, mesh_type, dimension, elements_neighbors_connection_);
            MeshFileHelpersFluent::nodeCoordinatesFluent( mesh_file, node_index, point_coordinates_3D_, text_line, dimension, point_coordinates);

            while (getline(mesh_file, text_line))
            {
                if (text_line.find("(13", 0) != string::npos && text_line.find(") (", 0) != string::npos)
                {
                 boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                 types_of_boundary_condition_.push_back(boundary_type);
                    while (getline(mesh_file, text_line))
                    {
                        if (text_line.find(")", 0) == string::npos)
                        {
                            Vecd nodes = MeshFileHelpers::nodeIndex(text_line);
                            Vec2d cells = MeshFileHelpers::cellIndex(text_line);
                            /*--- build up all topology---*/
                            bool check_neighbor_cell1 = 1;
                            bool check_neighbor_cell2 = 0;
                            for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                            {
                                if (mesh_type == 4)
                                {
                                    if (cells[check_neighbor_cell2] != 0)
                                    {
                                        MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);

                                        /*--- build up all connection data with element and neighbor and nodes---*/
                                        MeshFileHelpers::updateCellLists(cell_lists_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                        MeshFileHelpersFluent::updateBoundaryCellListsFluent(cell_lists_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                    }
                                    if (cells[check_neighbor_cell2] == 0)
                                    {
                                        break;
                                    }
                                }
                            }


                        }
                        else
                        break;
                    }
                }
                if (text_line.find("(12", 0) != string::npos && text_line.find("))", 0) != string::npos)
                break;
                if (text_line.find(")") != string::npos)
                continue;
            }
            cell_lists_.erase(cell_lists_.begin());
            
        }
        
    }
    //=================================================================================================//

        void readMeshFile_3d::getElementCenterCoordinates()
        {
            elements_center_coordinates_.resize(elements_nodes_connection_.size());
            elements_volumes_.resize(elements_nodes_connection_.size());
            for (std::vector<std::vector<long unsigned int>>::size_type element = 1; element != elements_nodes_connection_.size(); ++element)
            {
                Vecd center_coordinate = Vecd::Zero();
                MeshFileHelpers::cellCenterCoordinates(elements_nodes_connection_, element, point_coordinates_3D_, elements_center_coordinates_, center_coordinate);
                MeshFileHelpers::elementVolume(elements_nodes_connection_, element, point_coordinates_3D_, elements_volumes_);
            }
            elements_volumes_.erase(elements_volumes_.begin());
            elements_center_coordinates_.erase(elements_center_coordinates_.begin());
            elements_nodes_connection_.erase(elements_nodes_connection_.begin());
        }
        //=================================================================================================//

        void readMeshFile_3d::getMinimumDistanceBetweenNodes()
        {
            vector<Real> all_data_of_distance_between_nodes;
            all_data_of_distance_between_nodes.resize(0);
            MeshFileHelpers::faceArea(all_data_of_distance_between_nodes, elements_volumes_, cell_lists_, point_coordinates_3D_);
            auto min_distance_iter = std::min_element(all_data_of_distance_between_nodes.begin(), all_data_of_distance_between_nodes.end());
            if (min_distance_iter != all_data_of_distance_between_nodes.end())
            {
                min_distance_between_nodes_ = *min_distance_iter;
            }
            else
            {
                cout << "The array of all distance between nodes is empty " << endl;
            }
        }
        //=============================================================================================//

          //=================================================================================================//
          //=================================================================================================//
        ParticleGeneratorInFVM::ParticleGeneratorInFVM(SPHBody & sph_body, const StdLargeVec<Vecd> &positions, const StdLargeVec<Real> &elements_volumes)
            : ParticleGenerator(sph_body), elements_center_coordinates_(positions), elements_volumes_(elements_volumes) {}
        //=================================================================================================//
        void ParticleGeneratorInFVM::initializeGeometricVariables()
        {
            for (size_t particle_index = 0; particle_index != elements_center_coordinates_.size(); ++particle_index)
            {
                initializePositionAndVolumetricMeasure(elements_center_coordinates_[particle_index], elements_volumes_[particle_index]);
            }
        }
        //=================================================================================================//
        void BodyStatesRecordingInMeshToVtu::writeWithFileName(const std::string & sequence)
        {
            for (SPHBody* body : bodies_)
            {
                if (body->checkNewlyUpdated() && state_recording_)
                {
                    std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtu";
                    if (fs::exists(filefullpath))
                    {
                        fs::remove(filefullpath);
                    }
                    std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

                    MeshFileHelpers::vtuFileHeader(out_file);

                    MeshFileHelpers::vtuFileNodeCoordinates(out_file, nodes_coordinates_, elements_nodes_connection_);

                    MeshFileHelpers::vtuFileInformationKey(out_file);

                    MeshFileHelpers::vtuFileCellConnectivity(out_file, elements_nodes_connection_);

                    MeshFileHelpers::vtuFileOffsets(out_file, elements_nodes_connection_);

                    MeshFileHelpers::vtuFileTypeOfCell(out_file, elements_nodes_connection_);

                    //write Particle data to vtu file
                    out_file << "<CellData>\n";
                    body->writeParticlesToVtuFile(out_file);
                    out_file << "</CellData>\n";
                    // Write VTU file footer
                    out_file << "</Piece>\n";
                    out_file << "</UnstructuredGrid>\n";
                    out_file << "</VTKFile>\n";
                    out_file.close();
                }
                body->setNotNewlyUpdated();
            }
        }
    
}// namespace SPH