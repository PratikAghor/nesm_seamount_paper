* ```start_paths.m``` contains paths to required files, modify according to your system.
* ```save_vort_hslice.m``` saves horizontal slices of relative vorticity at a constant depth (-4500 m here, can be changed)
* ```save_vort_vslice.m``` saves vertical slices of relative vorticity at constant latitude or constant longitude
* ```save_uv_hslice_individually.m``` saves zonal (u) and meridional (v) components of velocity near the surface (at -100 m)
* ```read_plot_vort_hslice.m``` reads and plots horizontal slices of vorticity
* ```read_plot_vort_vslice_const_lat.m``` reads and plots vertical slices of vorticity at constant latitudes
* ```read_plot_avg_speed_hslice.py``` saves box-averaged timeseries of speed at (z = -100 m)
* ```aghor_hbbl``` directory contains scripts to plot vertical slices and horizontally-averaged timeseries of vertical diffusivity
* ```aghor_KmKe``` directory contains scripts to save and plot KmKe
* ```aghor_pv``` directory contains scripts to save and plot Ertel potential vorticity (PV)
* ```aghor_vebf``` directory contains scripts to save and plot VEBF
* Please read README.md files in the directories for instructions on running the postprocessing scripts
