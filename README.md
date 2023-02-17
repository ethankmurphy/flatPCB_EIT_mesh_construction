# flatPCB_EIT_mesh_construction
Codes to quickly construct FEM meshes from an STL output of a flat PCB design of electrodes

There are two main codes:
  * convert_CAD_elec_STL_to_elec_boundps_for_mesh.m
      * The user is prompted to pick the STL file
      * It automatically tries to find the flat direction
      * Then using ginput the user picks the electrode labeling
  * run_gmsh_mesh_from_CADebounds.m
      * Makes the outer surface triangulation (construct_outer_surf.m)
          * There’s an intermediate region – the size of it is picked automatically 
          * The shape of it is based on the shape of the electrodes and ensured to be at least a factor of 2x past the electrodes at any point
          * The meshing here is done using distmesh, https://github.com/ionhandshaker/distmesh
          * One needs just to download distmesh and add the path (lines 63-66 of run_gmsh_mesh_from_CADebounds.m)
      * Sets the h-values (local mesh density)
          * Assigned the h-electrode and h-subdomain based on factors of the minimum of the maximum width of all the electrodes 
          * Assigned the h-background based on a factor of the size of the subregion
       * Runs gmsh, then extracs the data, defines electrodes, and saves the data
          * gmsh is a free software, https://gmsh.info/
          * One needs to download/install the software and then supply the code the path to the executable (see lines 58-61 in run_gmsh_mesh_from_CADebounds.m)
  * check_plot_mesh.m: This last code is used to plot the mesh and qualitatively check everything is copacetic
  
There are two example STLs included in the folder, example_STLs,  
  * Electrode_Array_8ch_9o5mm v1.stl (see below with electrde labels)
![SMA_pic_with_elec_labels](https://user-images.githubusercontent.com/87721360/219694806-361511b4-9755-461c-b913-02c6684417d1.png)
  * A long electrode array, sarcopenia_US_EIT_Electrode_Array.stl
![sarcopenia_pic](https://user-images.githubusercontent.com/87721360/219695009-4daac293-d1b5-4175-a427-0af8d15cc288.png)

The output (check_plot_mesh.m) for the small electrode array (Electrode_Array_8ch_9o5mm v1) should look like this. 
![mesh_elec_bnds_Electrode_Array_8ch_9o5mm_v1_h0o1_0o2_5o2_3D_view](https://user-images.githubusercontent.com/87721360/219695693-b06294bc-19ec-4053-b8cd-797aa8b96627.png)

The h-values in run_gmsh_mesh_from_CADebounds.m automatically change based on the size of the PCB, but they can be manually changed. Additionally, one could adjust the maximum side or maximum depth limits of the resulting mesh if desired. 
