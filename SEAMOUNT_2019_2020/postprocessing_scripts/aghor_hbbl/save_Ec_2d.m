%-----------------------------------------
% save Ec(x, y) - depth integrated energy conversion rate from barotropic to baroclinic tides
% Ref. Pacific deep ocean circulation and ventilation controlled by tidal mixing awaay from the sea bottom
% Oka and Niwa (Nature Comm., 2013)
% Author: Pratik Aghor
%-----------------------------------------
% Ec(x, y) = depth integral of \bar(g rho' w_s) dz (eq. 7 in Oka & Niwa)
% rho' = rho - rho_0
% w_s = vertical velocity due to interaction of barotropic tidal flow and topography
% bar is time average
%-----------------------------------------
start_paths;
%--------------------------------------------------
indxRange = 469:3382;
[~, Nt] = size(indxRange);

%---------------------------------------------------
annual_mean_indxRange = 469:3382;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_file = strcat('seamount_2019_2020_rho_w_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

wb = ncread(annual_mean_file, 'wb');
rhob = ncread(annual_mean_file, 'rhob');

size(wb)
size(rhob)
%---------------------------------------------------
Ec_avg = zeros(M, L);

for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SEAMOUNT_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;

    zeta = pagetranspose(ncread(hisfile, 'zeta'));
    zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
    zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
    
    rhotmp = zeros(N, M, L);
        
    % construct 3d rho mat from get_rho.m
    for k=1:N
        vlevel = k;
        % tmp = get_rho(hisfile, gridfile, tindex, vlevel, coef);
        % size(tmp)
        [~, ~, ~, rhotmp(k, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    end
    size(rhotmp)

    ws = wb;
    rhoprime = (rhotmp - rhob);

    rpws = g.*(rhoprime.*ws);
    
    [Ec, h0] = vintegr2(abs(rpws), zw, zr, NaN, NaN); % perform vertical integration from bottom to top 
    Ec_avg = Ec_avg + (Ec);
end

Ec_avg = Ec_avg./Nt;

Ec_avg
%---------------------------------------------------
%---------------------------------------------------
%
% Save
%
% save data into a netcdf file
filename = strcat('seamount_2019_2020_Ec_2d_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
z_len = netcdf.defDim(ncid, 'Nz', N);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

vebf_varid = netcdf.defVar(ncid, 'Ec_2d', 'double', [y_len x_len]);
netcdf.putAtt(ncid, vebf_varid, 'description', 'time avg, depth integrated energy conversion rate from barotropic to baroclinic tides');
netcdf.putAtt(ncid, vebf_varid, 'units', 'kg s^-3');
netcdf.putAtt(ncid, vebf_varid, 'array dimensions', size(Ec_avg));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vebf_varid, Ec_avg);
% close netcdf file
netcdf.close(ncid);
disp('Done creating Ec_2d_avg!')
%--------------------------------------------------------
