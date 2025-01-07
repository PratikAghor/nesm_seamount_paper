%---------------------------------
% read and plot col avg EKE vs t
start_paths
%---------------------------------
indxRange = 469:3382; % entire year
vlevel = -4500;
[~, Nt] = size(indxRange);
% filename = strcat('nesm_2019_2020_kv_t_oka_niwa_nt_', ...
%     string(indxRange(1)), '_', string(indxRange(Nt)), '_vlevel_', string(vlevel), ...
%    '_old.nc');
filename = strcat('nesm_2019_2020_kv_t_oka_niwa_nt_', ...
    string(indxRange(1)), '_', string(indxRange(Nt)), '_vlevel_', string(vlevel), ...
   '_laurent_ferrari.nc');
k_bbl_t = ncread(filename, 'k_bbl_t');

filepath = strcat('../../../compare/k_bbl/kv_t_timeseries/');
k_bbl_filepath = strcat(filepath, 'nesm_2019_2020_kv_t_oka_niwa_nt_', ...
    string(indxRange(1)), '_', string(indxRange(Nt)),  '_vlevel_', string(vlevel), '.asc');
save(sprintf(k_bbl_filepath), 'k_bbl_t', '-ascii');
%--------------------------------------------------------------------------
% plot(eke_t, 'LineWidth', 2)