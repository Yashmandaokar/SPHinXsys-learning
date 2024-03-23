#include "MeshFileHelperFluent.h"
namespace SPH
{

    

    void MeshFileHelpersFluent::numberofNodesFluent(ifstream& mesh_file, size_t& number_of_points, string& text_line)
    {

        while (getline(mesh_file, text_line))

        {
            text_line.erase(0, 1);
            if (atoi(text_line.c_str()) == 10 && text_line.find("))", 0) != string::npos)
            {
            string text1(text_line);
            /*text_line.erase(3);
            string text2(text_line);
            text_line.erase(2);*/
            text1.erase(0, 8);
            Real last_position = text1.find_last_of(")");
            text1.erase(last_position - 3);
            number_of_points = stoi(text1, nullptr, 16);
            break;
            }
        }
    }

   void MeshFileHelpersFluent::numberofElementsFluent(ifstream& mesh_file, size_t& number_of_elements, string& text_line)
    {

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
                text1.erase(last_position - 3);
                number_of_elements = stoi(text1, nullptr, 16);
                break;
            }
        }
    }

    void MeshFileHelpersFluent::nodeCoordinatesFluent(ifstream& mesh_file, size_t& node_index, vector<vector<Real>> &point_coordinates_3D_, string& text_line, size_t& dimension, vector<vector<Real>>& point_coordinates)
    {
        while (getline(mesh_file, text_line))
        {
            text_line.erase(0, 1);
            if (atoi(text_line.c_str()) == 10 && text_line.find(") (", 0) != string::npos)
            {
                while (getline(mesh_file, text_line))
                { 
                    if (text_line.find("(", 0) == string::npos && text_line.find("))", 0) == string::npos && text_line.find(" ", 0) != string::npos)
                    {
                        if (node_index == 2019)
                        {
                            Real a = 1;
                        }
                        Real Coords[3] = { 0.0, 0.0, 0.0 };
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
                    text_line.erase(0, 1);
                    if (atoi(text_line.c_str()) == 11 || atoi(text_line.c_str()) == 13 || atoi(text_line.c_str()) == 12)
                    {
                        break;
                    }
                }
               
            }
        
            if ((atoi(text_line.c_str()) == 11 && text_line.find(") (") != string::npos) || (atoi(text_line.c_str()) == 13 && text_line.find(") (") != string::npos) || (atoi(text_line.c_str()) == 12 && text_line.find(") (") != string::npos))
            break;
        }
    }

    void MeshFileHelpersFluent::updateBoundaryCellListsFluent(vector<vector<vector<size_t>>>& cell_lists_, vector<vector<size_t>>& elements_nodes_connection_, Vecd nodes, Vec2d cells, bool& check_neighbor_cell1, bool& check_neighbor_cell2, size_t boundary_type)
    {

        if (cell_lists_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            cell_lists_[cells[check_neighbor_cell2]][1][0] = cells[check_neighbor_cell1];
            cell_lists_[cells[check_neighbor_cell2]][1][1] = boundary_type;
            cell_lists_[cells[check_neighbor_cell2]][1][2] = nodes[0];
            cell_lists_[cells[check_neighbor_cell2]][1][3] = nodes[1];
            cell_lists_[cells[check_neighbor_cell2]][1][4] = nodes[2];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
        if (cell_lists_[cells[check_neighbor_cell2]][0][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][1][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(0) && cell_lists_[cells[check_neighbor_cell2]][2][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1) && cell_lists_[cells[check_neighbor_cell2]][3][0] == static_cast<std::decay_t<decltype(elements_nodes_connection_[0][0])>>(-1))
        {
            cell_lists_[cells[check_neighbor_cell2]][2][0] = cells[check_neighbor_cell1];
            cell_lists_[cells[check_neighbor_cell2]][2][1] = boundary_type;
            cell_lists_[cells[check_neighbor_cell2]][2][2] = nodes[0];
            cell_lists_[cells[check_neighbor_cell2]][2][3] = nodes[1];
            cell_lists_[cells[check_neighbor_cell2]][2][4] = nodes[2];
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
        else
        {
            check_neighbor_cell1 = false;
            check_neighbor_cell2 = true;
            return;
        }
    }



}

