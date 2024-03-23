
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
            ofstream outfile("output.txt");
            if (!outfile) 
            {
                cerr << "Error: Unable to open file output.txt" << endl;
            }

            for (size_t i = 0; i < cell_lists_.size(); ++i) {
                for (size_t j = 0; j < cell_lists_[i].size(); ++j) {
                    for (size_t k = 0; k < cell_lists_[i][j].size(); ++k) {
                        outfile << cell_lists_[i][j][k] << " ";
                    }
                    outfile << endl;
                }
                outfile << endl;
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

    //=================================================================================================//
    BaseInnerRelationInFVM::BaseInnerRelationInFVM(RealBody &real_body, vector<vector<vector<size_t>>> data_inpute, vector<vector<Real>> nodes_coordinates)
        : BaseInnerRelation(real_body), real_body_(&real_body)
    {
        all_needed_data_from_mesh_file_ = data_inpute;
        nodes_coordinates_ = nodes_coordinates;
        subscribeToBody();
        resizeConfiguration();
    };
    //=================================================================================================//
    void BaseInnerRelationInFVM::resetNeighborhoodCurrentSize()
    {
        parallel_for(
            IndexRange(0, base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_),
            [&](const IndexRange &r)
            {
                for (size_t num = r.begin(); num != r.end(); ++num)
                {
                    inner_configuration_[num].current_size_ = 0;
                }
            },
            ap);
    }
    //=================================================================================================//
    void BaseInnerRelationInFVM::resizeConfiguration()
    {
        size_t updated_size = base_particles_.real_particles_bound_ + base_particles_.total_ghost_particles_;
        inner_configuration_.resize(updated_size, Neighborhood());
    }
            
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
    void NeighborBuilderInFVM::createRelation(Neighborhood &neighborhood, Real &distance,
                                            Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const
    {
        neighborhood.j_.push_back(j_index);
        neighborhood.r_ij_.push_back(distance);
        neighborhood.e_ij_.push_back(interface_normal_direction);
        neighborhood.dW_ijV_j_.push_back(dW_ijV_j);
        neighborhood.allocated_size_++;
    }
    //=================================================================================================//
    void NeighborBuilderInFVM::initializeRelation(Neighborhood &neighborhood, Real &distance,
                                                Real &dW_ijV_j, Vecd &interface_normal_direction, size_t j_index) const
    {
        size_t current_size = neighborhood.current_size_;
        neighborhood.j_[current_size] = j_index;
        neighborhood.dW_ijV_j_[current_size] = dW_ijV_j;
        neighborhood.r_ij_[current_size] = distance;
        neighborhood.e_ij_[current_size] = interface_normal_direction;
    }

    //=================================================================================================//
    InnerRelationInFVM::InnerRelationInFVM(RealBody &real_body, vector<vector<vector<size_t>>> data_inpute, vector<vector<Real>> nodes_coordinates)
        : BaseInnerRelationInFVM(real_body, data_inpute, nodes_coordinates), get_inner_neighbor_(&real_body){};
    //=================================================================================================//
    template <typename GetParticleIndex, typename GetNeighborRelation>
    void InnerRelationInFVM::searchNeighborsByParticles(size_t total_particles, BaseParticles &source_particles,
                                                        ParticleConfiguration &particle_configuration, GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation)
    {
        parallel_for(
            IndexRange(0, base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_),
            [&](const IndexRange &r)
            {
                StdLargeVec<Vecd> &pos_n = source_particles.pos_;
                StdLargeVec<Real> &Vol_n = source_particles.Vol_;
                for (size_t num = r.begin(); num != r.end(); ++num)
                {
                    size_t index_i = get_particle_index(num);
                    Vecd &particle_position = pos_n[index_i];
                    Real &Vol_i = Vol_n[index_i];

                    Neighborhood &neighborhood = particle_configuration[index_i];
                    for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != all_needed_data_from_mesh_file_[index_i].size(); ++neighbor)
                    {
                        size_t index_j = all_needed_data_from_mesh_file_[index_i][neighbor][0] - 1;
                        size_t boundary_type = all_needed_data_from_mesh_file_[index_i][neighbor][1];
                        size_t interface_node1_index = all_needed_data_from_mesh_file_[index_i][neighbor][2];
                        size_t interface_node2_index = all_needed_data_from_mesh_file_[index_i][neighbor][3];
                        size_t interface_node3_index = all_needed_data_from_mesh_file_[index_i][neighbor][4];
                        Vecd node1_position = Vecd(nodes_coordinates_[interface_node1_index][0], nodes_coordinates_[interface_node1_index][1], nodes_coordinates_[interface_node1_index][2]);
                        Vecd node2_position = Vecd(nodes_coordinates_[interface_node2_index][0], nodes_coordinates_[interface_node2_index][1], nodes_coordinates_[interface_node2_index][2]);
                        Vecd node3_position = Vecd(nodes_coordinates_[interface_node3_index][0], nodes_coordinates_[interface_node3_index][1], nodes_coordinates_[interface_node3_index][2]);
                        Vecd interface_area_vector1 = node1_position - node2_position;
                        Vecd interface_area_vector2 = node1_position - node3_position;
                        Vecd normal_vector = interface_area_vector1.cross(interface_area_vector2);
                        Real magnitude = normal_vector.norm();
                        Vecd normalized_normal_vector = normal_vector / magnitude;
                        Vecd particle_position = pos_n[index_i];
                        Vecd particle_position_j = pos_n[index_j];
                        Vecd node1_to_center_direction = particle_position - node1_position; 
                        if (node1_to_center_direction.dot(normalized_normal_vector) < 0)
                        {
                            normalized_normal_vector = -normalized_normal_vector; // vector pointing towards i
                        };
                        Real r_ij = 0; // we need r_ij to calculate the viscous force
                        // boundary_type == 2 means both of them are inside of fluid
                        if (boundary_type == 2)
                        {
                            r_ij = (particle_position - pos_n[index_j]).dot(normalized_normal_vector);
                        }
                        // boundary_type == 3 means fulid particle with wall boundary
                        if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 9) | (boundary_type == 10) | (boundary_type == 36) 
                            | (boundary_type == 5) | (boundary_type == 7))
                        {
                            r_ij = node1_to_center_direction.dot(normalized_normal_vector) * 2.0;
                        }
                        Real dW_ijV_j = (-0.5 * magnitude) / (2.0 * Vol_i);
                        get_neighbor_relation(neighborhood, r_ij, dW_ijV_j, normalized_normal_vector, index_j);
                    }
                }
            },
            ap);
    }
    //=================================================================================================//
    void InnerRelationInFVM::updateConfiguration()
    {
        resetNeighborhoodCurrentSize();
        searchNeighborsByParticles(base_particles_.total_real_particles_ + base_particles_.total_ghost_particles_,
                                base_particles_, inner_configuration_,
                                get_particle_index_, get_inner_neighbor_);
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