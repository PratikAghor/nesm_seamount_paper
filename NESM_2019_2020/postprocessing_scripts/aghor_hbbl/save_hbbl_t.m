%
% get density vs t and calculate height of the bottom boundary layer from rho
% Author: Pratik Aghor
%
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange = 469:3382; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);
%-------------------------------------------------------
rho_tmp = zeros(NumLayers, Ny, Nx);
% vert_rho = zeros(Nt, NumLayers);
hbbl_t = zeros(Nt, Ny, Nx);
drho_tol = 0.01; % to find hbbl, find where | rho- rho_bottom | < drho_tol, inside BBL
drho_truth_arr = zeros(NumLayers, 1); % 1 if | rho- rho_bottom | >= drho_tol, outside BBL

% plot options
z_depth = zeros(NumLayers, 1);
drho_vert_arr = zeros(NumLayers, 1); 
%% Calculate bottom boundary layer thickness

for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;

    % get density at each vlevel for this timestep nt
    for vlevel = 1:NumLayers
        [~, ~, ~, rho_tmp(vlevel, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    end
    
    size(rho_tmp)


    % find hbbl for each (i, j)
    for i = 1:Ny
        for j = 1:Nx
            z_depth(:) = squeeze(zr(:, i, j)); % depth 
            drho_vert_arr(:) =  rho_tmp(:, i, j) - rho_tmp(1, i, j); % difference from bottom-most density
        %%------------------------------------------
        % % Create figure
        % figure1 = figure(nt);
        % 
        % % Create axes
        % axes1 = axes('Parent',figure1,'YMinorTick','on',...
        %     'LineWidth',3,...
        %     'FontSize',24);
        % hold on
        % plot(drho_vert_arr, z_depth, 'LineWidth', 2);
        % 
        % xlabel('$\rho - \rho_0$','Interpreter','latex','FontSize',24); ylabel('$z$','Interpreter','latex','FontSize',24);
        % drho_vs_z  = strcat(drho_vs_z,'_nt_',num2str(nt));
        % exportgraphics(figure1, [drho_vs_z, '.pdf'], 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % close all
        %%------------------------------------------
        % find vlevel where drho_arr(vlevel) >= drho_tol from the bottom,
        % drho_truth_arr: inside BBL is 1, outside BBL is 0
        % also keep in mind zr goes from bottom to the surface
            for vlevel = 1:NumLayers
                if(abs(drho_vert_arr(vlevel)) < drho_tol)
                    drho_truth_arr(vlevel) = true; 
                end
            end

            [val, idx] = find(drho_truth_arr == 1);
            bbl_idx = max(idx) + 1; % this is the first vlevel outside of BBL
            
            % method 1:
            % hbbl = (interpolate between zr(i, j, bbl_idx) and zr(i, j, 1)) 
            % hbbl_t(nt, i, j) = 0.5*(zr(i, j, bbl_idx) + zr(i, j, bbl_idx)) - zr(i, j, 1); % height of BBL
            
            % method 2: linear interpolation: drho = a zr + b;
            % a = drho_vert_arr(bbl_idx)/(zr(i, j, bbl_idx) - zr(i, j, 1));
            % b = -a*zr(i, j, 1);
            % zr_bbl = (-drho_tol - b)/a; % at zr_bbl, drho = -drho_tol
            
            % method 2a: same as 2, using matlab function interp1
            x = [drho_vert_arr(bbl_idx) 0];
            v = [zr(bbl_idx, i, j) zr(1, i, j)];
            zr_bbl = interp1(x, v, -drho_tol);
            
            % method 3: just take the value at the next collocation point
            % zr_bbl = zr(i, j, bbl_idx);

            % calculate hbbl
            % hbbl_t(nt, i, j) = zr_bbl - zr(i, j, 1); % hbbl = zr(bbl) - zr(bottom)
            hbbl_t(nt, i, j) = zr_bbl + depth(i, j); % hbbl = zr(bbl) - zr(bottom)
            end
        end
    end
   %**********************************************************************
% for nt = 1:3
%     sprintf(strcat('loading', [dirHR, 'NESM1km_avg.%05d.nc'], 'file'), indxRange(1, nt))
%     fname = sprintf([dirHR, 'NESM1km_avg.%05d.nc'], indxRange(1, nt));
%     hisfile = fname;
% 
%     % get density at each vlevel for this timestep nt
%     for vlevel = 1:NumLayers
%         [~, ~, ~, rho_tmp(:, :, vlevel)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
%     end
% 
%     % find hbbl for each (i, j)
%     for i = 1:Ny
%         for j = 1:Nx
%             for k = 1:NumLayers
%                 if(abs(rho_tmp(i,j,k)-rho_tmp(i,j, 1)) > drho_tol)
%                     hbbl_t(nt, i, j) = - zr(i, j, k-1) + depth(i, j);
%                 end
%             end
%         end
%     end
% end
%**********************************************************************
%%
%% save hbbl
% save data into a netcdf file
hbbl_t_filename = strcat('nesm_2019_2020_hbbl_t_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');

ncid = netcdf.create(hbbl_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(hbbl_t_filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

hbbl_t_varid = netcdf.defVar(ncid, 'hbbl_t', 'double', [t_len y_len x_len]);
netcdf.putAtt(ncid, hbbl_t_varid, 'description', 'time series of hbbl for lon-lat values');
netcdf.putAtt(ncid, hbbl_t_varid, 'units', 'meter');
netcdf.putAtt(ncid, hbbl_t_varid, 'array dimensions', size(hbbl_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
netcdf.putVar(ncid, hbbl_t_varid, hbbl_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating hbbl_t for lon-lat values!')

%%------------------------------------------

