* ```save_rho_w_annual_mean.m```: saves annual mean of rho and w, or if it has already been saved to calculate VEBF, read annual mean from ```../aghor_vebf/``` directory.
* ```save_Ec_2d.m``` saves depth integrated energy conversion rate from barotropic to baroclinic tides according to Oka and Niwa (Nature Comm. 2013) 
* ```save_hbbl_t.m``` saves the height of the bottom boundary layer (BBL) for each lat-lon location.
* ```save_sigma4_t_bbl_slices.m``` saves potential density at 4000 dbar using TEOS10 (can change to 2000 dbar as well using ```gsw_sigma2``` instead of ```gsw_sigma4```)
* ```save_kv_t_vslice_laurent_ferrari.m``` saves vertical slices of vertical diffusivity as a modified expression of (St Laurent et. al., 2002), according to Oka and Niwa (Nature Comm. 2013) and (Ferrari et. al., 2016)
* ```save_kv_t_hslice_avg_laurent_ferrari.m``` saves horizontal slice of vertical diffusivity at a constant depth of -4500 m. 
