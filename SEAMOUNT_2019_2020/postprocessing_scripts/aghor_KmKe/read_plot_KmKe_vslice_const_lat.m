%-------------------------------------------------------
% plot averaged rho vlsice over given timeframe (typically 7 days)
% ref: Winter Mixed Layer Restratification Induced by Vertical Eddy
% Buoyancy Flux in the Labrador Sea, Li et al., GRL
% vertical eddy buoyancy flux (vebf) = <w'b'>

% Goal: to obtain b'w', with b = -g*(rho - rho0)/rho0; 
% primes denote deviations from a time average

% plot b'w'(z, t) 
% plot <b'w'>_z -> z-avgeraged time series 
%%
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
% case_type ='instant'; % 'instant' or 'avg'
case_type ='avg'; % 'instant' or 'avg'

if(strcmp(case_type, 'instant'))
    nt = 3040;
elseif(strcmp(case_type, 'avg'))
    % indxRange = 469:3382; % entire year
    % indxRange = 712:768; % April 01 - April 07, 2019
    % indxRange = 1040:1096; % May 12 - May 19, 2019
    % indxRange = 2104:2160; % Sept 22 - Sept 29, 2019
    indxRange = 3016:3072; % Jan 14 - Jan 20, 2020
    nt0=indxRange(1);
    [~, Nt] = size(indxRange);
end
%-------------------------------------------------------
% only to get the dates in the title!
vlevel = -4500;
vort_indxRange = 0:3382; 
[~, vort_Nt] = size(vort_indxRange);
nt0=vort_indxRange(1);

vort_t_filename = strcat('../seamount_2019_2020_vort_t_hslice_vlevel_', string(vlevel(1)), '_nt_',string(vort_indxRange(1)), '_', ...
    string(vort_indxRange(vort_Nt)), '.nc');
t_arr = ncread(vort_t_filename, 't_arr');
t_i = t_arr(indxRange(1) - nt0);
t_e = t_arr(indxRange(Nt) - nt0);
date_i = char(string(datetime(2017, 12, 31, 23, 29, t_i)));
date_e = char(string(datetime(2017, 12, 31, 23, 29, t_e)));

date_i = date_i(1:6);
date_e = date_e(1:11);

date = strcat(date_i, '-', date_e);
%-------------------------------------------------------

%-------------------------------------------------------
%%
lon_idx_arr = [120, 140, 165, 180];
X = (repmat(lon_rho_vec, NumLayers, 1));
%%
%-------------------------------------------------------
if(strcmp(case_type, 'instant'))
    filename = strcat('seamount_2019_2020_instant_KmKe_vslice_const_lat_nt_',string(nt), '.nc');

elseif(strcmp(case_type, 'avg'))
    filename = strcat('seamount_2019_2020_KmKe_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
        string(indxRange(Nt)), '.nc');
end

hrs_vslice_arr = (ncread(filename, 'hrs_vslice'));
vrs_vslice_arr = (ncread(filename, 'vrs_vslice'));
kmke_vslice_arr = (ncread(filename, 'kmke_vslice'));

lat_arr = (ncread(filename, 'lat'));
lat_idx_arr = (ncread(filename, 'lat_idx'));
%-------------------------------------------------------
%%
% mkdir vort_plots;
for i = 1:length(lat_idx_arr)
    lat_idx = (lat_idx_arr(i)); % find(abs(lat_rho_vec-lat_arr(1))<1e-3); % find idx of lat in lat_rho_vec
    Z = (squeeze(zr(:,lat_idx,:)));
    %-------------------------------------
    %% plotting HRS vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    hrs_vslice = squeeze(hrs_vslice_arr(:, i, :));
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, hrs_vslice); 
    colormap(cmocean('haline')); colorbar; 
    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
    
    for j = 1:length(lon_idx_arr)
            xline(squeeze(X(:, lon_idx_arr(j))));
    end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])

    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    if(strcmp(case_type, 'instant'))
        figname  = [plots_path, 'rs_plots/hrs/seamount_2019_2020_instant_hrs_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(nt));
    elseif(strcmp(case_type, 'avg'))
        figname  = [plots_path, 'rs_plots/hrs/seamount_2019_2020_hrs_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
        string(indxRange(Nt)));
    end

    % exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %
    % save for all i
    % if(i == 1)
        % for j = 1:length(lon_idx_arr)
        %     lon_idx = lon_idx_arr(j);
        %     % save file to compare/vebf
        %     data = [squeeze(Z(:, lon_idx)), squeeze(hrs_vslice(:, lon_idx))];
        %     filename = '../../../compare/KmKe/vslice/KmKe/seamount_2019_2020_KmKe';
        %     filename = strcat(filename, '_lat_', string(lat_arr(i)));
        %     filename = strcat(filename, '_lon_', string(X(1, lon_idx)));
        %     filename = strcat(filename, '_nt_', string(indxRange(1)), '_', ...
        %     string(indxRange(Nt)));
        %     filename = strcat(filename, '.asc');
        %     dlmwrite(filename, data, 'delimiter', '\t');
        % end
    % end
    %-------------------------------------
    %% plotting VRS vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    vrs_vslice = squeeze(vrs_vslice_arr(:, i, :));
    % KmKe_vslice = hrs_vslice + vrs_vslice;
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, vrs_vslice); 
    colormap(cmocean('haline')); colorbar; 
    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
    
    % for j = 1:length(lon_idx_arr)
    %         xline(squeeze(X(:, lon_idx_arr(j))));
    % end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])

    title(ax1, string(date));

    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    if(strcmp(case_type, 'instant'))
        figname  = [plots_path, 'rs_plots/vrs/seamount_2019_2020_instant_vrs_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(nt));
    elseif(strcmp(case_type, 'avg'))
        figname  = [plots_path, 'rs_plots/vrs/seamount_2019_2020_vrs_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
        string(indxRange(Nt)));
    end
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %-------------------------------------
    %% plotting KmKe vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    kmke_vslice = squeeze(kmke_vslice_arr(:, i, :));
    % KmKe_vslice = hrs_vslice + vrs_vslice;
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, kmke_vslice); 
    % colormap(cmocean('haline')); colorbar;
    colormap('parula'); cb = colorbar;
    cb.FontSize = 20;
    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
     
    % for j = 1:length(lon_idx_arr)
    %         xline(squeeze(X(:, lon_idx_arr(j))));
    % end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])
    title(ax1, string(date));
    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    if(strcmp(case_type, 'instant'))
        figname  = [plots_path, 'rs_plots/KmKe/seamount_2019_2020_instant_KmKe_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(nt));
    elseif(strcmp(case_type, 'avg'))
        figname  = [plots_path, 'rs_plots/KmKe/seamount_2019_2020_KmKe_vslice'];
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
        string(indxRange(Nt)));
    end
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %-------------------------------------
end
%%
