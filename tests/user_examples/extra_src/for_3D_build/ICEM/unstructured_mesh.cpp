
#include "unstructured_mesh.h"
namespace SPH
{
    //=================================================================================================//
    void readMeshFile_3d::getDataFromMeshFile3d()
    {
        ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
        mesh_file.open(full_path_);
        if (mesh_file.fail())
        {
            cout << "Error:Check if the file exists." << endl;
        }
        string text_line;
        /*--- Read the dimension of the problem ---*/
        size_t dimension(0);
        while (getline(mesh_file, text_line))
        {
            text_line.erase(0, 1);
            text_line.erase(0, 2);
            istringstream value(text_line);
            if (text_line.find("3", 0) != string::npos)
            {
                dimension = atoi(text_line.c_str());
                break;
            }
        }
        /*--- Read the node data (index is starting from zero) ---*/
        size_t number_of_points(0);
        while (getline(mesh_file, text_line))
        {
            text_line.erase(0, 1);
            string text1(text_line);
            text_line.erase(3);
            string text2(text_line);
            text_line.erase(2);
            if (atoi(text_line.c_str()) == 10 && text1.find("))", 0) != string::npos)
            {
                text1.erase(0, 8);
                Real last_position = text1.find_last_of(")");
                text1.erase(last_position - 5);
                number_of_points = stoi(text1, nullptr, 16);
                break;
            }
        };
        point_coordinates_3D_.resize(number_of_points);
        for (std::vector<std::vector<double>>::size_type point = 0; point != point_coordinates_3D_.size(); ++point)
        {
            point_coordinates_3D_[point].resize(dimension);
        }

        /* Prepare our data structure for the point coordinates. */
        point_coordinates.resize(dimension);
        for (std::size_t k = 0; k < dimension; k++)
        {
            point_coordinates[k].reserve(number_of_points);
        }
        size_t node_index = 0;
        while (getline(mesh_file, text_line))
        {
            if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos)
            {
                Real Coords[3] = { 0.0, 0.0, 0.0 };
                if (text_line.find(" ", 0) != string::npos)
                {
                    if (dimension == 2)
                    {
                        size_t devide_position = text_line.find_first_of(" ");
                        string x_part = text_line;
                        string y_part = text_line;
                        string x_coordinate_string = x_part.erase(devide_position);
                        string y_coordinate_string = y_part.erase(0, devide_position);

                        istringstream streamx, streamy;
                        streamx.str(x_coordinate_string);
                        streamy.str(y_coordinate_string);
                        streamx >> Coords[0];
                        streamy >> Coords[1];
                        point_coordinates_3D_[node_index][0] = Coords[0];
                        point_coordinates_3D_[node_index][1] = Coords[1];
                        ++node_index;
                    }
                    if (dimension == 3)
                    {
                        size_t first_devide_position = text_line.find_first_of(" ");
                        size_t last_devide_position = text_line.find_last_of(" ");
                        string x_part = text_line;
                        string y_part = text_line;
                        string z_part = text_line;
                        string x_coordinate_string = x_part.erase(first_devide_position);
                        string y_coordinate_string = y_part.erase(last_devide_position);
                        y_coordinate_string = y_coordinate_string.erase(0, first_devide_position);
                        string z_coordinate_string = z_part.erase(0, last_devide_position);
                        istringstream streamx, streamy, streamz;
                        streamx.str(x_coordinate_string);
                        streamy.str(y_coordinate_string);
                        streamz.str(z_coordinate_string);
                        streamx >> Coords[0];
                        streamy >> Coords[1];
                        streamz >> Coords[2];
                        point_coordinates_3D_[node_index][0] = Coords[0];
                        point_coordinates_3D_[node_index][1] = Coords[1];
                        point_coordinates_3D_[node_index][2] = Coords[2];
                        ++node_index;
                    }
                    for (std::size_t iDim = 0; iDim != dimension; ++iDim)
                    {
                        point_coordinates[iDim].push_back(Coords[iDim]);
                    }
                }
            }
            if (text_line.find("))", 0) != string::npos)
            {
                break;
            }
        }

        /*--- Read the elements of the problem (13 (id f1 f2 type 0) (---*/
        /** differnet boundary conditions
         * bc-type==2, interior boundary condition.
         * bc-type==3, wall boundary condition.
         * bc-type==4, Pressure Inlet boundary condition
         * bc-type==5, Pressure Outlet boundary condition
         * bc-type==7, Symmetry boundary condition
         * bc-type==9, pressure-far-field boundary condition.
         * bc-type==a, Velocity Inlet boundary condition.
         * bc-type==c, Periodic boundary condition.
         * bc-type==e, porous jumps boundary condition.
         * bc-type==14,Mass Flow Inlet boundary condition.
         * Note that Cell0 means boundary condition.
         * mesh_type==4, unstructured mesh.
         * mesh_type==6, structured mesh.
         * cell_lists_
         * {[(neighbor_cell_index, bc_type, node1_of_face, node2_of_face, node3_of_face), (....), (.....)], []..... }.
         * {inner_neighbor1, inner_neighbor2,inner_neighbor3, inner_neighbor4 ..... }.
         */
        size_t boundary_type(0);
        size_t number_of_elements(0);
        size_t mesh_type = 4;
        while (getline(mesh_file, text_line))
        {
            text_line.erase(0, 1);
            string text1(text_line);
            text_line.erase(3);
            string text2(text_line);
            text_line.erase(2);
            if (atoi(text_line.c_str()) == 12)
            {
                text1.erase(0, 8);
                Real last_position = text1.find_last_of(")");
                text1.erase(last_position - 5);
                number_of_elements = stoi(text1, nullptr, 16);
                break;
            }
        };
        cell_lists_.resize(number_of_elements + 1);
        for (std::size_t a = 0; a != number_of_elements + 1; ++a)
        {
            cell_lists_[a].resize(mesh_type);
            for (std::vector<std::vector<long unsigned int>>::size_type b = 0; b != cell_lists_[a].size(); ++b)
            {
                cell_lists_[a][b].resize(dimension + 2);
                for (std::vector<long unsigned int>::size_type c = 0; c != cell_lists_[a][b].size(); ++c)
                {
                    cell_lists_[a][b][c] = -1;
                }
            }
        }
        /*--- reinitialize the number of elements ---*/
        elements_nodes_connection_.resize(number_of_elements + 1);
        elements_neighbors_connection_.resize(number_of_elements + 1);
        cell_lists_.resize(number_of_elements + 1);
        for (std::size_t element = 0; element != number_of_elements + 1; ++element)
        {
            elements_nodes_connection_[element].resize(4);
            for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
            {
                elements_nodes_connection_[element][node] = -1;
            }
        }

        /*--- find the elements lines ---*/
        while (getline(mesh_file, text_line))
        {
            if (text_line.find("(13", 0) != string::npos && text_line.find(")(", 0) != string::npos)
            {
                /*--- find the type of boundary condition ---*/
                size_t position = text_line.find(")", 0);
                text_line = text_line.erase(0, position - 4);
                text_line = text_line.erase(2);
                boundary_type = stoi(text_line, nullptr, 16);
                types_of_boundary_condition_.push_back(boundary_type);
                while (getline(mesh_file, text_line))
                {
                    if (text_line.find(")", 0) == string::npos)
                    {
                        /*--- find the node1 between two cells ---*/
                        string node1_string_copy = text_line;
                        size_t first_devide_position = text_line.find_first_of(" ", 0);
                        string node1_index_string = node1_string_copy.erase(first_devide_position);
                        size_t node1_index_decimal = stoi(node1_index_string, nullptr, 16) - 1;

                        /*--- find the node2 between two cells---*/
                        string node2_string = text_line;
                        string node3_string = text_line;
                        node2_string = node2_string.erase(0, first_devide_position + 1);
                        node3_string = node2_string;
                        size_t second_devide_position = node2_string.find_first_of(" ", 0);
                        node2_string.erase(second_devide_position);
                        size_t node2_index_decimal = stoi(node2_string, nullptr, 16) - 1;

                        /*--- find the node3 between two cells---*/

                        node3_string = node3_string.erase(0, second_devide_position + 1);
                        size_t third_devide_position = node3_string.find_first_of(" ", 0);
                        node3_string.erase(third_devide_position);
                        size_t node3_index_decimal = stoi(node3_string, nullptr, 16) - 1;

                        Vecd nodes = Vecd(node1_index_decimal, node2_index_decimal, node3_index_decimal);

                        /*--- find the cell1---*/
                        string cell1_string = text_line;
                        cell1_string = cell1_string.erase(0, first_devide_position + 1);
                        cell1_string = cell1_string.erase(0, second_devide_position + 1);
                        cell1_string = cell1_string.erase(0, third_devide_position + 1);
                        size_t fourth_devide_position = cell1_string.find_first_of(" ", 0);
                        cell1_string.erase(fourth_devide_position);
                        size_t cell1_index_decimal = stoi(cell1_string, nullptr, 16);

                        /*--- find the cell2---*/
                        string cell2_string = text_line;
                        cell2_string = cell2_string.erase(0, first_devide_position + 1);
                        cell2_string = cell2_string.erase(0, second_devide_position + 1);
                        cell2_string = cell2_string.erase(0, third_devide_position + 1);
                        cell2_string.erase(0, fourth_devide_position + 1);
                        size_t cell2_index_decimal = stoi(cell2_string, nullptr, 16);
                        Vec2d cells = Vec2d(cell1_index_decimal, cell2_index_decimal);
                        if (cell2_index_decimal == 1091)
                        {
                           Real y = 1;
                        }
                       
                        /*--- build up all topology---*/
                        bool check_neighbor_cell1 = 1;
                        bool check_neighbor_cell2 = 0;
                        for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                        {
                            while (true)
                            {
                                if (mesh_type == 4)
                                {
                                    if (cells[check_neighbor_cell2] != 0)
                                    {
                                        /*--- build up connection with element and nodes only---*/
                                        for (int node = 0; node != nodes.size(); ++node)
                                        {
                                            if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][1] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][2] != nodes[node] && elements_nodes_connection_[cells[check_neighbor_cell2]][3] != nodes[node])
                                            {
                                                if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && (elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1)) && elements_nodes_connection_[cells[check_neighbor_cell2]][2] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][3] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                                {
                                                    elements_nodes_connection_[cells[check_neighbor_cell2]][3] = nodes[node];
                                                }
                                                if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][1] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][2] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                                {
                                                    elements_nodes_connection_[cells[check_neighbor_cell2]][2] = nodes[node];
                                                }
                                                if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && elements_nodes_connection_[cells[check_neighbor_cell2]][1] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                                {
                                                    elements_nodes_connection_[cells[check_neighbor_cell2]][1] = nodes[node];
                                                }
                                                if (elements_nodes_connection_[cells[check_neighbor_cell2]][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                                {
                                                    elements_nodes_connection_[cells[check_neighbor_cell2]][0] = nodes[node];
                                                }

                                            }
                                            else
                                                continue;
                                        }
                                        /*--- build up all connection data with element and neighbor and nodes---*/
                                        if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][2][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][3][0] != cells[check_neighbor_cell1])
                                        {
                                            if (cell_lists_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                /*--- inner neighbor index---*/
                                                cell_lists_[cells[check_neighbor_cell2]][0][0] = cells[check_neighbor_cell1];
                                                /*--- boundary type---*/
                                                cell_lists_[cells[check_neighbor_cell2]][0][1] = boundary_type;
                                                /*--- nodes of a face---*/
                                                cell_lists_[cells[check_neighbor_cell2]][0][2] = nodes[0];
                                                cell_lists_[cells[check_neighbor_cell2]][0][3] = nodes[1];
                                                cell_lists_[cells[check_neighbor_cell2]][0][4] = nodes[2];

                                                check_neighbor_cell1 = false;
                                                check_neighbor_cell2 = true;
                                                break;
                                            }
                                            if (cell_lists_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                cell_lists_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
                                                cell_lists_[cells[check_neighbor_cell2]][1][1] = boundary_type;
                                                cell_lists_[cells[check_neighbor_cell2]][1][2] = nodes[0];
                                                cell_lists_[cells[check_neighbor_cell2]][1][3] = nodes[1];
                                                cell_lists_[cells[check_neighbor_cell2]][1][4] = nodes[2];
                                                check_neighbor_cell1 = false;
                                                check_neighbor_cell2 = true;
                                                break;
                                            }
                                            if (cell_lists_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                                                cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                                                cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                                                cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                                                cell_lists_[cells[check_neighbor_cell2]][2][4] = nodes[2];
                                                check_neighbor_cell1 = false;
                                                check_neighbor_cell2 = true;
                                                break;
                                            }
                                            if (cell_lists_[cells[check_neighbor_cell2]][0][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][1][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][2][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                            {
                                                cell_lists_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
                                                cell_lists_[cells[check_neighbor_cell2]][3][1] = boundary_type;
                                                cell_lists_[cells[check_neighbor_cell2]][3][2] = nodes[0];
                                                cell_lists_[cells[check_neighbor_cell2]][3][3] = nodes[1];
                                                cell_lists_[cells[check_neighbor_cell2]][3][4] = nodes[2];
                                                check_neighbor_cell1 = false;
                                                check_neighbor_cell2 = true;
                                                break;
                                            }
                                        }
                                    
                                        if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][1][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            cell_lists_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
                                            cell_lists_[cells[check_neighbor_cell2]][3][1] = boundary_type;
                                            cell_lists_[cells[check_neighbor_cell2]][3][2] = nodes[0];
                                            cell_lists_[cells[check_neighbor_cell2]][3][3] = nodes[1];
                                            cell_lists_[cells[check_neighbor_cell2]][3][4] = nodes[2];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
                                            cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
                                            cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
                                            cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
                                            cell_lists_[cells[check_neighbor_cell2]][2][4] = nodes[2];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        } 
                                        if (cell_lists_[cells[check_neighbor_cell2]][0][0] != cells[check_neighbor_cell1] && cell_lists_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][2][0] != static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
                                        {
                                            cell_lists_[cells[check_neighbor_cell2]][3][0] = cells[check_neighbor_cell1];
                                            cell_lists_[cells[check_neighbor_cell2]][3][1] = boundary_type;
                                            cell_lists_[cells[check_neighbor_cell2]][3][2] = nodes[0];
                                            cell_lists_[cells[check_neighbor_cell2]][3][3] = nodes[1];
                                            cell_lists_[cells[check_neighbor_cell2]][3][4] = nodes[2];
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                        
                                        else
                                        {
                                            check_neighbor_cell1 = false;
                                            check_neighbor_cell2 = true;
                                            break;
                                        }
                                    }
                                    if (cells[check_neighbor_cell2] == 0)
                                    {
                                        break;
                                    }
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
    //=================================================================================================//

    void readMeshFile_3d::getElementCenterCoordinates()
    {
        elements_center_coordinates_.resize(elements_nodes_connection_.size());
        elements_volumes_.resize(elements_nodes_connection_.size());
        Real Total_volume = 0;
        for (std::vector<std::vector<long unsigned int>>::size_type element = 1; element != elements_nodes_connection_.size(); ++element)
        {
            Vecd center_coordinate = Vecd::Zero();
            for (std::vector<long unsigned int>::size_type node = 0; node != elements_nodes_connection_[element].size(); ++node)
            {
                center_coordinate += Vecd(point_coordinates_3D_[elements_nodes_connection_[element][node]][0] / 4.0,
                    point_coordinates_3D_[elements_nodes_connection_[element][node]][1] / 4.0, point_coordinates_3D_[elements_nodes_connection_[element][node]][2] / 4.0);
            }
            elements_center_coordinates_[element] = center_coordinate;

            // get nodes position
            using Vec4d = Eigen::Matrix<Real, 4, 1>;
            Vec4d nodes = Vec4d(elements_nodes_connection_[element][0], elements_nodes_connection_[element][1], elements_nodes_connection_[element][2], elements_nodes_connection_[element][3]);
            Vecd node1_coordinate = Vecd(point_coordinates_3D_[nodes[0]][0], point_coordinates_3D_[nodes[0]][1], point_coordinates_3D_[nodes[0]][2]);
            Vecd node2_coordinate = Vecd(point_coordinates_3D_[nodes[1]][0], point_coordinates_3D_[nodes[1]][1], point_coordinates_3D_[nodes[1]][2]);
            Vecd node3_coordinate = Vecd(point_coordinates_3D_[nodes[2]][0], point_coordinates_3D_[nodes[2]][1], point_coordinates_3D_[nodes[2]][2]);
            Vecd node4_coordinate = Vecd(point_coordinates_3D_[nodes[3]][0], point_coordinates_3D_[nodes[3]][1], point_coordinates_3D_[nodes[3]][2]);
            // get each line length
         /*   Real first_side_length = (node1_coordinate - node2_coordinate).norm();
            Real second_side_length = (node1_coordinate - node3_coordinate).norm();
            Real third_side_length = (node1_coordinate - node4_coordinate).norm();
            Real fourth_side_length = (node2_coordinate - node3_coordinate).norm();
            Real fifth_side_length = (node2_coordinate - node4_coordinate).norm();
            Real sixth_side_length = (node3_coordinate - node4_coordinate).norm();*/
            // half perimeter
            //Real half_perimeter = (first_side_length + second_side_length + third_side_length + fourth_side_length + fifth_side_length + sixth_side_length) / 2.0;     

            // Create the 3x3 matrix M
            Mat3d M;
            M << node2_coordinate - node1_coordinate, node3_coordinate - node1_coordinate, node4_coordinate - node1_coordinate;
            Real determinant = abs(M.determinant());
            Real element_volume = (static_cast<double>(1) / 6) * determinant;
            if (element == 53151)
            {
               Real j = 1;
            }
            Total_volume += element_volume;
            //Real element_volume =
            //   pow(half_perimeter * (half_perimeter - first_side_length) * (half_perimeter - second_side_length) * (half_perimeter - third_side_length), 0.5);
            elements_volumes_[element] = element_volume;
            
        }
        std::ofstream outfile("elements_nodes_connection.txt");
        if (outfile.is_open()) 
        {
            // Iterate through the vector and write each element to the file
            for (const auto& element : elements_nodes_connection_)
            {
                for (const int nod : element)
                {
                    outfile << nod << " ";
                }
                outfile << "\n"; // Newline after each element
            }
            outfile.close();
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
    for (size_t element_index = 0; element_index != elements_volumes_.size(); ++element_index)
    {
        for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != cell_lists_[element_index].size(); ++neighbor)
        {
            /*if (element_index == 36150)
            {
              Real  a = 1;
            }*/
            size_t interface_node1_index = cell_lists_[element_index][neighbor][2];
            size_t interface_node2_index = cell_lists_[element_index][neighbor][3];
            size_t interface_node3_index = cell_lists_[element_index][neighbor][4];
            Vecd node1_position = Vecd(point_coordinates_3D_[interface_node1_index][0], point_coordinates_3D_[interface_node1_index][1], point_coordinates_3D_[interface_node1_index][2]);
            Vecd node2_position = Vecd(point_coordinates_3D_[interface_node2_index][0], point_coordinates_3D_[interface_node2_index][1], point_coordinates_3D_[interface_node2_index][2]);
            Vecd node3_position = Vecd(point_coordinates_3D_[interface_node3_index][0], point_coordinates_3D_[interface_node3_index][1], point_coordinates_3D_[interface_node3_index][2]);
            Vecd interface_area_vector1 = node2_position - node1_position;
            Vecd interface_area_vector2 = node3_position - node1_position;
            Vecd area_vector = interface_area_vector1.cross(interface_area_vector2);
            Real triangle_area = 0.5*area_vector.norm();
            //Real interface_area_size = interface_area_vector1.norm();
            all_data_of_distance_between_nodes.push_back(triangle_area);
        }
    }
    std::ofstream outfile("Cell_list.txt");
    for (size_t i = 0; i < cell_lists_.size(); ++i) 
    {
        for (size_t j = 0; j < cell_lists_[i].size(); ++j) 
        {
            for (size_t k = 0; k < cell_lists_[i][j].size(); ++k) 
            {
                outfile << cell_lists_[i][j][k] << " ";
            }
            outfile << std::endl;  // Start a new line for each set of 5 elements
        }
    }
    outfile.close();
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
    ParticleGeneratorInFVM::ParticleGeneratorInFVM(SPHBody &sph_body, const StdLargeVec<Vecd> &positions, const StdLargeVec<Real> &elements_volumes)
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
void BodyStatesRecordingInMeshToVtu::writeWithFileName(const std::string &sequence)
{
    for (SPHBody *body : bodies_)
    {
        if (body->checkNewlyUpdated() && state_recording_)
        {
            // TODO: we can short the file name by without using SPHBody
            std::string filefullpath = io_environment_.output_folder_ + "/SPHBody_" + body->getName() + "_" + sequence + ".vtu";
            if (fs::exists(filefullpath))
            {
                fs::remove(filefullpath);
            }
            std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
            // begin of the XML file
            //out_file << "<?xml version=\"1.0\"?>\n";
            out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
            out_file << "<UnstructuredGrid>\n";
            out_file << "<FieldData>\n";
            out_file << "<DataArray type=\"Int32\" Name=\"ispatch\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";
            out_file << "0\n";
            out_file << "</DataArray>\n";
            out_file << "</FieldData>\n";
             
            // Write point data
            out_file << "<Piece NumberOfPoints=\"" << nodes_coordinates_.size() << "\" NumberOfCells=\"" << elements_nodes_connection_.size() << "\">\n";
            out_file << "<PointData>\n";
            out_file << "</PointData>\n";
            out_file << "<CellData>\n";
            out_file << "</CellData>\n";
            out_file << "<Points>\n";
            out_file << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"4.124318125496386\">\n";

            size_t total_nodes = nodes_coordinates_.size();
            for (size_t node = 0; node != total_nodes; ++node)
            {
                out_file << nodes_coordinates_[node][0] << " " << nodes_coordinates_[node][1] << " " << nodes_coordinates_[node][2] << " \n";
                /*out_file << nodes_coordinates_[node][0] << " " << nodes_coordinates_[node][1] << " 0.0\n";*/
            }

            out_file << "<InformationKey name=\"L2_NORM_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
            out_file << "<Value index=\"0\">\n";
            out_file << "0\n";
            out_file << "</Value>\n";
            out_file << "<Value index=\"1\">\n";
            out_file << "4.1243181255\n";
            out_file << "</Value>\n";
            out_file << "</InformationKey>\n";
            out_file << "<InformationKey name=\"L2_NORM_FINITE_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
            out_file << "<Value index=\"0\">\n";
            out_file << "0\n";
            out_file << "</Value>\n";
            out_file << "<Value index=\"1\">\n";
            out_file << "4.1243181255\n";
            out_file << "</Value>\n";
            out_file << "</InformationKey>\n";
            out_file << "</DataArray>\n";
            out_file << "</Points>\n";
            
            
            out_file << "<Cells>\n";
            // Write Cell data
            out_file << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"11138\">\n";
/*
            for (int i = 0; i < elements_nodes_connection_.size; ++i) 
            {
            for (int j = 0; j < 4; ++j) 
            {
                outfile << enc[i+1][j] << " ";
            }
                std::cout << std::endl;
            }*/

            for (const auto& cell : elements_nodes_connection_)
            {
                for (const auto& vertex : cell)
                {
                    out_file << vertex << " ";
                }
                out_file << "\n";
            }

            out_file << "</DataArray>\n";
            out_file << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"212604\">\n";

            size_t offset = 0;
            for (const auto& face : elements_nodes_connection_)
            {
                offset += face.size();
                out_file << offset << " ";
            }


            size_t type = 10;
            out_file << "\n</DataArray>\n";

            //Specifies type of cell 10 = tetrahedral
            out_file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"10\" RangeMax=\"10\">\n";
           for (const auto& types : elements_nodes_connection_)
            {
                for (const auto& vertex : types)
                {
                    out_file << type << " ";
                }
                out_file << "\n";
            }
            // Write face attribute data
            out_file << "</DataArray>\n";
            out_file << "</Cells>\n";
            //write Particle data to vtu file
            out_file << "<CellData>\n";
            body->writeParticlesToVtuFile(out_file);
            out_file << "</CellData>\n";
            // Write file footer
            out_file << "</Piece>\n";
            out_file << "</UnstructuredGrid>\n";
            out_file << "</VTKFile>\n";
            
            out_file.close();
        }
        body->setNotNewlyUpdated();
    }
}


}// namespace SPH