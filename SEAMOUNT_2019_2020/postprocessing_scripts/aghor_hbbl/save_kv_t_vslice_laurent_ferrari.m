%----------------------------------------------------------
% save vertical slices of avg kv according to Ferrari et. al. + Oka-Niwa
% Ref: Turning ocean upside down - Ferrari et. al. (eq. 14)
% modify to: kv = -vebf/N^2 + Gamma*q*Ec(x, y)*F(z)/(rho*N^2) 
% Pratik Aghor
%----------------------------------------------------------
%--------------------------------------------------
start_paths;
%--------------------------------------------------
% indxRange = 469:470;
indxRange = 3016:3072; % Jan 14, Jan 20, 2020
vlevel = -4500;
dz = 50;
vlevel2 = vlevel - dz;

[~, Nt] = size(indxRange);
h_decay = 500; % 500 m, Oka & Niwa, St Laurent, etc.
q = 0.3334;
Gamma = 0.2;
Kb = 1e-5; % background diffusivity
%---------------------------------------------------
%---------------------------------------------------
vebf_file = strcat('../aghor_vebf/seamount_2019_2020_VEBF_3d_nt_',string(indxRange(1)), '_', string(indxRange(Nt)),'.nc');
vebf_3d = ncread(vebf_file, 'VEBF3D');
% zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
% zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
[vebf_2d, h0] = vintegr2((vebf_3d), zw, zr, NaN, NaN);
%---------------------------------------------------
% hbbl
hbbl_indxRange = 469:3382
hbbl_nt0 = hbbl_indxRange(1); % 469
[~, hbbl_Nt] = size(hbbl_indxRange);
hbbl_filename =  strcat('seamount_2019_2020_hbbl_t_nt_',string(hbbl_indxRange(1)), '_', string(hbbl_indxRange(hbbl_Nt)),'.nc');
hbbl_t = ncread(hbbl_filename, 'hbbl_t');
%---------------------------------------------------
Ec_indxRange = 469:3382;
% Ec_indxRange=3016:3072;
[~, Ec_Nt] = size(Ec_indxRange);
Ec_file = strcat('seamount_2019_2020_Ec_2d_nt_',string(Ec_indxRange(1)), '_', string(Ec_indxRange(Ec_Nt)),'.nc');
Ec_2d = ncread(Ec_file, 'Ec_2d');
Ec_2d = abs(Ec_2d);
%---------------------------------------------------
dz_r = zr(2:end,:,:)-zr(1:end-1,:,:);
time_arr = zeros(Nt, 1);
%---------------------------------------------------
lat_idx_arr = ([90, 112, 135]);
lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
kv_t_vslice_const_lat = zeros(Nt, N, length(lat_arr), L);
%---------------------------------------------------
lon_idx_arr = ([150, 158, 166, 172]);
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j)); % const lon vals to save vslice at
end
kv_t_vslice_const_lon = zeros(Nt, N, Ny, length(lon_arr));
%---------------------------------------------------
for nt = 1:Nt
    
    hbbl_nt = indxRange(nt) + nt - hbbl_nt0;
    hbbl = squeeze(hbbl_t(hbbl_nt, :, :));

    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;

    time = ncread(fname, 'time');
    time_arr(nt) = time;
    
    w = pagetranspose(ncread(hisfile, 'w'));
    w = shiftdim(w, 2);

    rhotmp = zeros(N, Ny, Nx);
    % construct 3d rho mat from get_rho.m
    for k=1:N
        vlevel = k;
        [~, ~, ~, rhotmp(k, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    end
    % get 3d temp, salt in (Nz, Ny, Nx) form
    % temp = pagetranspose(ncread(hisfile, 'temp'));
    % salt = pagetranspose(ncread(hisfile, 'salt'));

    % temp = shiftdim(temp, 2);
    % salt = shiftdim(salt, 2);
    
    % rhotmp  = gsw_sigma4(salt, temp); % get potential density
    % rhotmp = rhotmp + rho0;

    drho_dz = zeros(N, M, L);

    drho_dz(2:end, :, :) = squeeze((rhotmp(2:N, :, :) - rhotmp(1:N-1, :, :)))./...
        (dz_r);
    drho_dz(1, :, :) = drho_dz(2, :, :);

    N2_tmp = (-g./rho0).*drho_dz;
    
    b_zz = zeros(N, M, L);
    b_zz(2:end, :, :) =  squeeze((N2_tmp(2:N, :, :) - N2_tmp(1:N-1, :, :)))./...
        (dz_r);
    b_zz(1, :, :) = b_zz(2, :, :);
    % calculate kv_3d_t
    kv_3d_tmp = zeros(N, M, L);
    for k=1:N
        for i=1:Ny
            for j =1:Nx
                z = zr(k, i, j);
                zb_minus_z = depth(i, j) + z;
                h_decay = hbbl(i, j); % testing non constant h_decay
                F_z_at_vlevel = exp(- (zb_minus_z)/h_decay)./(h_decay*(1 - exp(-depth(i, j)/h_decay)));
		
	    %    h_decay_vebf_term = depth(i, j);
		%    F_z_at_vlevel_vebf_term = exp(- (zb_minus_z)/h_decay_vebf_term)./((1 - exp(-depth(i, j)/h_decay_vebf_term)));
                	
                % kv_3d_tmp(k, i, j) = 0.5*(Kb + ...
		%			+ (Gamma*(q*(Ec_2d(i, j)/rhotmp(k, i, j)))*F_z_at_vlevel/(N2_tmp(k, i, j))  ) ...
		%			- ( vebf_3d(k, i, j)*F_z_at_vlevel_vebf_term/(N2_tmp(k, i, j))  ));
		
		if(vebf_3d(k, i, j) >= 0)
			k0 = Kb;
		else
			k0 = -vebf_3d(k, i, j)./N2_tmp(k, i, j);
		end

		kv_3d_tmp(k, i, j) = k0 + Gamma*(q*(Ec_2d(i, j)/rhotmp(k, i, j)))*F_z_at_vlevel/(N2_tmp(k, i, j));

		if(kv_3d_tmp(k, i, j) <= 0 )
		   kv_3d_tmp(k, i, j) = Kb;
		end
            end
        end
    end

    % kv_3d_tmp = Kb + ((Gamma).*vebf_3d)./N2_tmp;

    % now take const_lat_vslices
    for i = 1:length(lat_arr)
	lat_val = lat_arr(i);
        lat_idx = lat_idx_arr(i);
	kv_t_vslice_const_lat(nt, :, i, :) = squeeze(kv_3d_tmp(:, lat_idx, :));
    end
    
    % save const_lon_vslices
    for j = 1:length(lon_arr)
	lon_val = lon_arr(j);
	lon_idx = lon_idx_arr(j);
	kv_t_vslice_const_lon(nt, :, :, j) = squeeze(kv_3d_tmp(:, :, lon_idx));
    end
end
%---------------------------------------------------
% save
% const lat vslice
% save data into a netcdf file
filename = strcat('seamount_2019_2020_kv_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'_laurent_ferrari.nc');
ncid = netcdf.create(filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
y_len = netcdf.defDim(ncid, 'Nlat', length(lat_arr));
x_len = netcdf.defDim(ncid, 'Nx', Nx);
depth_len = netcdf.defDim(ncid, 'Nz', N);

netcdf.close(ncid);
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

lat_idx_varid = netcdf.defVar(ncid, 'lat_idx', 'double', [y_len]);
netcdf.putAtt(ncid, lat_idx_varid, 'description', 'indices in lat_rho_vec for const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(y_len));

lat_varid = netcdf.defVar(ncid, 'lat', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical kv slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(y_len));

vort_varid = netcdf.defVar(ncid, 'kv_t_vslice', 'double', [t_len depth_len y_len x_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of vertical slices of kv_t at fixed lat values');
netcdf.putAtt(ncid, vort_varid, 'units', 'm^2/s');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(kv_t_vslice_const_lat));

% close define mode
netcdf.endDef(ncid);

% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, vort_varid, kv_t_vslice_const_lat);
% close netcdf file
netcdf.close(ncid);
%---------------------------------------------------
%% const lon vort_vslice
% save data into a netcdf file
filename = strcat('seamount_2019_2020_kv_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'_laurent_ferrari.nc');
ncid2 = netcdf.create(filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid2, 'Nt', Nt);
y_len = netcdf.defDim(ncid2, 'Ny', Ny);
x_len = netcdf.defDim(ncid2, 'Nlon', length(lon_arr));
depth_len = netcdf.defDim(ncid2, 'Nz', N);

netcdf.close(ncid2);
% define variables and attributes
ncid2 = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid2);

time_varid = netcdf.defVar(ncid2, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid2, time_varid, 'description', 'time array');
netcdf.putAtt(ncid2, time_varid, 'units', 's');


lon_idx_varid = netcdf.defVar(ncid2, 'lon_idx', 'double', [x_len]);
netcdf.putAtt(ncid2, lon_idx_varid, 'description', 'indices in lon_rho_vec for const lon values where vertical epv slice is taken');
netcdf.putAtt(ncid2, lon_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lon_idx_varid, 'array dimensions', size(x_len));

lon_varid = netcdf.defVar(ncid2, 'lon', 'double', [x_len]);
netcdf.putAtt(ncid2, lon_varid, 'description', 'const lon values where vertical kv slice is taken');
netcdf.putAtt(ncid2, lon_varid, 'units', '--');
netcdf.putAtt(ncid, lon_varid, 'array dimensions', size(x_len));

vort_varid = netcdf.defVar(ncid2, 'kv_t_vslice', 'double', [t_len depth_len y_len x_len]);
netcdf.putAtt(ncid2, vort_varid, 'description', 'time series of normalized vertical slices of kv at fixed lon values');
netcdf.putAtt(ncid2, vort_varid, 'units', 'm^2/s');
netcdf.putAtt(ncid2, vort_varid, 'array dimensions', size(kv_t_vslice_const_lon));

% close define mode
netcdf.endDef(ncid2);

% put values
netcdf.putVar(ncid2, time_varid, time_arr); % start from 0
netcdf.putVar(ncid2, lon_idx_varid, lon_idx_arr);
netcdf.putVar(ncid2, lon_varid, lon_arr);
netcdf.putVar(ncid2, vort_varid, kv_t_vslice_const_lon);
% close netcdf file
netcdf.close(ncid2);
%%
%---------------------------------------------------
disp('Done creating kv_t vslices according to Ferrari et. al.!')
%---------------------------------------------------
