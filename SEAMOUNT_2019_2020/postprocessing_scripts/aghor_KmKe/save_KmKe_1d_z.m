% 
% save KmKe vs z for horizontal averages over given lan-lon boxes
% Pratik Aghor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Assume mean and KmKe_3d_zlevs has been calculated and saved.

clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
% indxRange = 1992:2800; % What time indices do you need?
% indxRange = 469:3382; % entire year
% indxRange = 712:768; % Apr 01 - Apr 07, 2019
% indxRange = 1040:1096; % May 12 - May 19, 2019
% indxRange = 2080:2144; % Sept 2019
% indxRange = 2104:2160; % Sept 22 - Sept 29, 2019

indxRange = 3016:3072; % Jan 14 - Jan 20, 2020
% indxRange = 3080:3144;

nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);

% KmKe_3d_indxRange = indxRange; % 469:3382;
% [~, Nt] = size(KmKe_3d_indxRange);
KmKe_3d_file = strcat('seamount_2019_2020_KmKe_3d_zlevs_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

zlevs = ncread(KmKe_3d_file, 'zlevs');
KmKe3D_zlevs = ncread(KmKe_3d_file, 'KmKe3D_zlevs');
HRS3D_zlevs = ncread(KmKe_3d_file, 'HRS3D_zlevs');
VRS3D_zlevs = ncread(KmKe_3d_file, 'VRS3D_zlevs');

KmKe1D = zeros(N, 2);
HRS1D = zeros(N, 2);
VRS1D = zeros(N, 2);
%--------------------------------------------------------------------------
% define box idx to take horizontal averages
% box_arr = ["N", "NE", "W", "SW", "S", "SE", "E", "NW"]
box_arr = ["big_N", "big_S"]
for i = 1:length(box_arr)
    box = string(box_arr(i))
    [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);

    %--------------------------------------------------------------------------
    for k = 1:N
        KmKe1D(k, 2) = mean(squeeze(KmKe3D_zlevs(k, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max)), "all", "omitnan");
        HRS1D(k, 2) = mean(squeeze(HRS3D_zlevs(k, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max)), "all", "omitnan");
        VRS1D(k, 2) = mean(squeeze(VRS3D_zlevs(k, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max)), "all", "omitnan");
    end
    
    KmKe1D(:, 1) = zlevs;
    HRS1D(:, 1) = zlevs;
    VRS1D(:, 1) = zlevs;
    % replace Nan values with 0
    % KmKe1D(isnan(KmKe1D))=0; 
    % HRS1D(isnan(HRS1D))=0;
    % VRS1D(isnan(VRS1D))=0;
    %--------------------------------------------------------------------------
    filepath = strcat('../../../compare/KmKe/vslice/');
    KmKe_filepath = strcat(filepath, 'KmKe/seamount_2019_2020_kmke_z_', box, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
    HRS_filepath = strcat(filepath, 'hrs/seamount_2019_2020_hrs_z_' , box, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
    VRS_filepath = strcat(filepath, 'vrs/seamount_2019_2020_vrs_z_', box, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
    save(sprintf(HRS_filepath), 'HRS1D', '-ascii');
    save(sprintf(VRS_filepath), 'VRS1D', '-ascii');
    save(sprintf(KmKe_filepath), 'KmKe1D', '-ascii');
    %--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
% vertical line profiles at const lat-lon vs z
lat_idx_arr = ([90, 112, 135]);
lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
%--------------------------------------------------------------------------
lon_idx_arr = ([150, 158, 166, 172]);
% lon_idx_arr = ([172]);
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j)); % const lon values to save vslice at
end
latsec = lat_rho_vec;
%--------------------------------------------------------
%--------------------------------------------------------------------------
for i = 1:length(lat_arr)
    lat_idx = lat_idx_arr(i);
    for j = 1:length(lon_arr)
        lon_idx = lon_idx_arr(j);
        for k = 1:N
            KmKe1D(k, 2) = (KmKe3D_zlevs(k, lat_idx, lon_idx));
            HRS1D(k, 2) = (HRS3D_zlevs(k, lat_idx, lon_idx));
            VRS1D(k, 2) = ((VRS3D_zlevs(k, lat_idx, lon_idx)));
            KmKe1D(:, 1) = zlevs;
            HRS1D(:, 1) = zlevs;
            VRS1D(:, 1) = zlevs;
            
            % replace Nan values with 0
            % KmKe1D(isnan(KmKe1D))=0; 
            % HRS1D(isnan(HRS1D))=0;
            % VRS1D(isnan(VRS1D))=0;

            filepath = strcat('../../../compare/KmKe/vslice/');
            KmKe_filepath = strcat(filepath, 'KmKe/seamount_2019_2020_kmke_z_lat_', string(lat_arr(i)), '_lon_', string(lon_arr(j)), '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
            HRS_filepath = strcat(filepath, 'hrs/seamount_2019_2020_hrs_z_lat_' , string(lat_arr(i)), '_lon_', string(lon_arr(j)), '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
            VRS_filepath = strcat(filepath, 'vrs/seamount_2019_2020_vrs_z_lat_', string(lat_arr(i)), '_lon_', string(lon_arr(j)), '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)),  '.asc');
            save(sprintf(HRS_filepath), 'HRS1D', '-ascii');
            save(sprintf(VRS_filepath), 'VRS1D', '-ascii');
            save(sprintf(KmKe_filepath), 'KmKe1D', '-ascii');
            %--------------------------------------------------------------

        end
    end
end

%--------------------------------------------------------------------------
