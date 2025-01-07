%--------------------------------------------------
% save box averaged kv timeseries at a constant depth
% using Oka-Niwa and St. Laurent formulation 
% Ref. Pacific deep ocean circulation and ventilation controlled by tidal mixing awaay from the sea bottom
% Oka and Niwa (Nature Comm., 2013)
% Author: Pratik Aghor

% assume Ec_2d has been saved
%--------------------------------------------------
start_paths;
%--------------------------------------------------
indxRange = 469:3382;
vlevel = -4500;
dz = 50; 
vlevel2 = vlevel - dz;

[~, Nt] = size(indxRange);
h_decay = 100; % 500 m, Oka & Niwa, St Laurent, etc.
q = 0.3334;
Gamma = 0.2;
Kb = 1e-5; % background diffusivity
%---------------------------------------------------
vebf_file = strcat('../aghor_vebf/seamount_2019_2020_VEBF_3d_nt_',string(indxRange(1)), '_', string(indxRange(Nt)),'.nc');
vebf_3d = ncread(vebf_file, 'VEBF3D');
% zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
% zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
% [vebf_2d, h0] = vintegr2((vebf_3d), zw, zr, NaN, NaN);
%---------------------------------------------------
%---------------------------------------------------
% hbbl
hbbl_indxRange = 469:3382
hbbl_nt0 = hbbl_indxRange(1); % 469
[~, hbbl_Nt] = size(hbbl_indxRange);
hbbl_filename =  strcat('seamount_2019_2020_hbbl_t_nt_',string(hbbl_indxRange(1)), '_', string(hbbl_indxRange(hbbl_Nt)),'.nc');
hbbl_t = ncread(hbbl_filename, 'hbbl_t');
%---------------------------------------------------

Ec_indxRange = 469:3382;
[~, Ec_Nt] = size(Ec_indxRange);
Ec_file = strcat('seamount_2019_2020_Ec_2d_nt_',string(Ec_indxRange(1)), '_', string(Ec_indxRange(Ec_Nt)),'.nc');
Ec_2d = ncread(Ec_file, 'Ec_2d');
%---------------------------------------------------
k_2d = zeros(Ny, Nx);
k_t_box_avg = zeros(Nt, 1);
for nt = 1:Nt

    hbbl_nt = nt + indxRange(1) - hbbl_nt0
    hbbl = squeeze(hbbl_t(hbbl_nt, :, :));

    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;

    % w = pagetranspose(ncread(hisfile, 'w'));

    % w = shiftdim(w, 2);

       
    % construct 3d rho mat from get_rho.m
    [~, ~, ~, rho_hslice] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    % [~, ~, ~, N2_hslice_from_bvf_eos] = get_bvf(hisfile, gridfile, tindex, vlevel, coef); % already N2 (1/s^2), see bvf_eos.m, 
    [~, ~, ~, rho_hslice_2] = get_rho(hisfile, gridfile, tindex, vlevel2, coef);
    N2_hslice = (rho_hslice - rho_hslice_2)./(dz);
    N2_hslice = (-g/rho0).*N2_hslice;
    
    % get vebf_hslice
    vebf_hslice = vinterp(vebf_3d, zr, vlevel);

    for i = 1:Ny
	for j = 1:Nx
        zb_minus_z = depth(i, j) + vlevel;
	h_decay = hbbl(i, j); % non constant h_decay
	F_z_at_vlevel = exp(- (zb_minus_z)/h_decay)./(h_decay*(1 - exp(-depth(i, j)/h_decay)));
	
	if(N2_hslice(i, j) < 1e-7)
		N2_hslice(i, j) = 5e-5; % to avoid unrealistic values of Kv, see Oka & Niwa (2013)
	end
        

   	% k_2d(i, j) = -(Gamma*vebf_hslice(i, j)/N2_hslice(i, j)) + (Gamma*q*Ec_2d(i, j)*F_z_at_vlevel/(rho_hslice(i, j)*N2_hslice(i, j)) );

	if(vebf_hslice(i, j) >= 0)
        	k0 = Kb;
        else
        	k0 = -vebf_hslice(i, j)./N2_hslice(i, j);
        end

        k_2d(i, j) = k0 + (Gamma*q*Ec_2d(i, j)*F_z_at_vlevel/(rho_hslice(i, j)*N2_hslice(i, j)) );

        if(k_2d(i, j) <= 0 )
                   k_2d(i, j) = Kb;
        end

	end
    end
    
    k_t_box = k_2d(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);
    k_t_box(isinf(k_t_box)|isnan(k_t_box)) = 0;
    disp('flag')
    size(k_t_box)

    % k_t_box_new = k_t_box(isfinite(k_t_box));
    k_t_box_avg(nt) = mean(k_t_box, 'all', 'omitnan'); % mean(mean(k_t_box, 'omitnan'), 'omitnan');
    % set a cap to avoid unrealistic values
    % cap_val = 4e-3;
    % if(k_t_box_avg(nt) >= cap_val)
    %    k_t_box_avg(nt) = cap_val;
    % end
    sprintf('kv_t_box_avg = %.7e', k_t_box_avg(nt) ) 
end
%%------------------------------------------
% save timeseries of box avg into a netcdf file
filename = strcat('seamount_2019_2020_kv_t_oka_niwa_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '_vlevel_', string(vlevel(1)), '_laurent_ferrari.nc')

ncid = netcdf.create(filename,'CLOBBER');
%% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% z_len = netcdf.defDim(ncid, 'Nz', Nz);
netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

k_t_varid = netcdf.defVar(ncid, 'k_bbl_t', 'double', [t_len]);
netcdf.putAtt(ncid, k_t_varid, 'description', 'timeseries of box averaged near bottom diffusivity');
netcdf.putAtt(ncid, k_t_varid, 'units', 'm^2s^-1');
netcdf.putAtt(ncid, k_t_varid, 'array dimensions', size(k_t_box_avg));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, k_t_varid, k_t_box_avg);
% close netcdf file
netcdf.close(ncid);
disp('Done creating near bottom diffusivity!')
%%------------------------------------------

