%-------------------------------------------------------
% read and plot 2d vertical slices of instantaneous vorticity
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
nt=3045;
%-------------------------------------------------------
%-------------------------------------------------------
%%
% lon_idx_arr = [120, 140, 165, 180];
X = (repmat(lon_rho_vec, NumLayers, 1));
%%
%-------------------------------------------------------
filename = strcat('seamount_2019_2020_instant_vort_vslice_const_lat_nt_',string(nt),'.nc');

vort_vslice_arr = (ncread(filename, 'vort_t_vslice'));

lat_arr = (ncread(filename, 'lat'));
lat_idx_arr = (ncread(filename, 'lat_idx'));
%-------------------------------------------------------


%%
% mkdir vort_plots;
for i = 1:length(lat_idx_arr)
    lat_idx = (lat_idx_arr(i)); % find(abs(lat_rho_vec-lat_arr(1))<1e-3); % find idx of lat in lat_rho_vec
    Z = (squeeze(zr(:,lat_idx,:)));
    %-------------------------------------
    %% plotting vort vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    vort_vslice = squeeze(vort_vslice_arr(:, i, :));
    % KmKe_vslice = hrs_vslice + vrs_vslice;
    zMin = -1; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, vort_vslice); 
    colormap(cmocean('balance')); colorbar; 
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

    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    figname  = [plots_path, 'vort_plots/vslice/const_lat/seamount_2019_2020_vort_vslice'];
    figname = strcat(figname, '_lat_', string(lat_arr(i)));
    figname = strcat(figname, '_nt_', string(nt));
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %-------------------------------------
end
%%
