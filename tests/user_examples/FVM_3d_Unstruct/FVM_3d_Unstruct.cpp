/**
 * @file 	Yash_3d_FVM_double_mach_reflection_Unstruct.cpp
 * @brief 	This is the compressible test for the realization of FVM in the SPHinXsys.
 * @details We consider a double mach reflection case.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "FVM_3d_Unstruct.h"
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
    // read data from ANSYS mesh.file
    readMeshFile_3d read_mesh_data(double_mach_reflection_Unstruct_mesh_fullpath);
    //----------------------------------------------------------------------
   //	Build up the environment of a SPHSystem.
   //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody wave_block(sph_system, makeShared<WaveBody>("WaveBody"));
    wave_block.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_another, heat_capacity_ratio);
    wave_block.generateParticles<ParticleGeneratorInFVM>(read_mesh_data.elements_center_coordinates_, read_mesh_data.elements_volumes_);
    wave_block.addBodyStateForRecording<Real>("Density");
    wave_block.addBodyStateForRecording<Real>("Pressure");
    /** Initial condition and register variables*/
    SimpleDynamics<DMFInitialCondition> initial_condition(wave_block);
    GhostCreationFromMesh ghost_creation(wave_block, read_mesh_data.cell_lists_, read_mesh_data.point_coordinates_3D_);

    // Visualization in FVM with date in cell.
    BodyStatesRecordingInMeshToVtu write_real_body_states(
        io_environment, sph_system.real_bodies_, read_mesh_data.elements_nodes_connection_, read_mesh_data.point_coordinates_3D_);
    write_real_body_states.writeToFile(0);
    return 0;
}