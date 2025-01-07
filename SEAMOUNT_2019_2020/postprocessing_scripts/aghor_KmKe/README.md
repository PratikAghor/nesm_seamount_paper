* Ref: Topographic and Mixed Layer Submesoscale Currents in the Near-Surface Southwestern Tropical Pacific, Sreenivasan et. al, JPO, 2017
* Modified from Diagnostic_tools/KmKe

#### for hslice
* first save_KmKe_3d: saves HRS, VRS, KmKe 3d matrices of time averages
* save_KmKe 2d takes horizontal slices at a given vlevel for HRS, VRS and KmKe
* read_plot_KmKe_2d plots hslices of HRS, VRS, KmKe at the given vlevel
#### for vslice
* first save_KmKe_3d: saves HRS, VRS, KmKe 3d matrices of time averages
* then save_KmKe_vslice_const_lat
* then read_plot_KmKe_vslice_const_lat 

#### convert to zlevs
* save_KmKe_3d.m saves HRS, VRS, KmKe 3d matrices of time averages in sigma coordinates
* To convert these 3d matrices to 3d zlevs, do save_KmKe_3d_zlevs.m
* Then take horizontal averages and save z variations using save_KmKe_1d_z.m 
