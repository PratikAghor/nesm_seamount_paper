* Ref: Topographic and Mixed Layer Submesoscale Currents in the Near-Surface Southwestern Tropical Pacific, Sreenivasan et. al, JPO, 2017
* Modified from Diagnostic_tools/KmKe

#### for hslice
* first save_vebf_3d: saves VEBF 3d matrices of time averages
* save_vebf 2d takes horizontal slices at a given vlevel
* read_plot_vebf_2d_hslice plots hslices of VEBF at the given vlevel
#### for vslice
* first save_vebf_3d: saves VEBF 3d matrices of time averages
* then save_vebf_vslice_const_lat
* then read_plot_vebf_vslice_const_lat 

#### convert to zlevs
* ```save_vebf_3d.m``` saves VEBF 3d matrices of time averages in sigma coordinates
* To convert these 3d matrices to 3d zlevs, do ```save_vebf_3d_zlevs.m```
* Then take horizontal averages and save z variations using ```save_vebf_1d_z.m``` 
