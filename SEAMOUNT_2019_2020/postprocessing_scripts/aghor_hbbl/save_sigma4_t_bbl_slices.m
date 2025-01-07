%--------------------------------------------------------------------
% save sigma4 potential density at 4000 dbar pressure using TEOS 10
% save 3 layers of rho_pot corresponding to sigma levels 1, 2, 3 (3 bottom most layers)
% need only 3 layers to calculate x and z derivatives in the expression for diapycnal diffusivity at BBL 
% Author: Pratik Aghor
%--------------------------------------------------------------------
start_paths;
%--------------------------------------------------------------------
indxRange=469:3382;
[~, Nt] = size(indxRange);


Nz = 3; % save only 3 sigma levels (bottom most)
sigma4_t_3d = zeros(Nt, Nz, Ny, Nx);

for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    
    % get 3d temp, salt in (Nz, Ny, Nx) form
    temp = pagetranspose(ncread(hisfile, 'temp'));
    salt = pagetranspose(ncread(hisfile, 'salt'));
    
    temp = shiftdim(temp, 2);
    salt = shiftdim(salt, 2);

    % get only the bottom-most 3 layers
    temp = temp(1:Nz, :, :);
    salt = salt(1:Nz, :, :);

    % size(temp)
    % size(salt)
    
    sigma4_t_3d(nt, :, :, :)  =  gsw_sigma4(salt, temp);  
end
%-----------------------------------------
% save data into a netcdf file
hbbl_t_filename = strcat('seamount_2019_2020_sigma4_t_3d_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');

ncid = netcdf.create(hbbl_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
z_len = netcdf.defDim(ncid, 'Nz', Nz);
netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(hbbl_t_filename,'WRITE');
netcdf.reDef(ncid);

sigma4_t_varid = netcdf.defVar(ncid, 'sigma4_t_3d', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, sigma4_t_varid, 'description', 'sigma4_t 3d matrix in terrain following sigma--coordinates');
netcdf.putAtt(ncid, sigma4_t_varid, 'units', 'kgm^-3');
netcdf.putAtt(ncid, sigma4_t_varid, 'array dimensions', size(sigma4_t_3d));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, sigma4_t_varid, sigma4_t_3d);
% close netcdf file
netcdf.close(ncid);
disp('Done creating sigma4_t_3d in terrain--following sigma coordinates!')
%%------------------------------------------
