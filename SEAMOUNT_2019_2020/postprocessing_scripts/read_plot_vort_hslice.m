%-------------------------------------------------------
%%
start_paths

%-------------------------------------------------------
%%
%%Loop to calculate averages if need (modify as needed)
% indxRange = 0:1690; % What time indices do you need?

indxRange = 0:3382; 
uv_indxRange=0:3382;


vlevel = ([-4500]);
nt0=indxRange(1);
uv_nt0=uv_indxRange(1);

% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
[~, uv_Nt] = size(uv_indxRange);

ke_arr = zeros(Nt, NumLayers); % get ke at each time at different depths
time_arr = zeros(Nt, 1);
%--------------------------------------------------------
tskip = 1;
vlevel_skip = 50;
vort_t_filename = strcat('seamount_2019_2020_vort_t_hslice_vlevel_', string(vlevel(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
u_t_filename = strcat('seamount_2019_2020_u_t_hslice_vlevel_-100_nt_',string(uv_indxRange(1)), '_', ...
    string(uv_indxRange(uv_Nt)),'.nc');

v_t_filename = strcat('seamount_2019_2020_v_t_hslice_vlevel_-100_nt_',string(uv_indxRange(1)), '_', ...
    string(uv_indxRange(uv_Nt)),'.nc');


% 
t_arr = ncread(vort_t_filename, 't_arr');
% vlevel = ncread(vort_t_filename, 'vlevel');
vort_t = ncread(vort_t_filename, 'vort_t');
u_t = ncread(u_t_filename, 'u_t');
v_t = ncread(v_t_filename, 'v_t');

f0 = 0.909*1e-4; % Coriolis freq 
%%
%-------------------------------------------------------

%-------------------------------------------------------
% mkdir vort_plots;
%% plotting
for nt = 2130:tskip:2133
    sprintf('for loop indx = %d', nt)
    vort_t_nt = indxRange(1) - nt0 + nt;
    sprintf('vort_t arr indx = %d', vort_t_nt)
    sprintf('fig indx = %d', vort_t_nt + nt0 - 1)
    for k = 1:1
        %%% plotting
        figure1 = figure(nt);
        [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        % ax1 = axes('LineWidth', 3, 'FontSize', 24);

        m_proj('miller','long', lonlim,'lat', latlim);
        % h1 = m_pcolor(lon_rho, lat_rho, log10(R));
       

        hold on;
        vort = squeeze(vort_t(vort_t_nt, :, :, k))/f0;
        u_hslice = squeeze(u_t(vort_t_nt, :, :));
        v_hslice = squeeze(v_t(vort_t_nt, :, :));
        zMin = -0.5; % min(min(vort)); 
        zMax = 0.5; % max(max(vort));
        h1 = m_image(lon_rho_vec, lat_rho_vec, vort);
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        % colormap(ax1,b2r(zMin,zMax));  colorbar;
        % colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
        % colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
        % freezeColors; hold on;
       
        xyskip = 15;

        m_quiver(lon_rho(1:xyskip:end, 1:xyskip:end), lat_rho(1:xyskip:end, 1:xyskip:end), ...
         u_hslice(1:xyskip:end, 1:xyskip:end), v_hslice(1:xyskip:end, 1:xyskip:end), ...
         'color',[0 0 0]);
       
        hold on;
        % add a reference arrow using m_vec
        [hpv5, htv5] = m_vec(100, -62, 40.25, 20, 0, 'k', 'key', '0.2 m/s');
        % [hpv5, htv5] = m_vec(100, -61.5, 37.4, 20, 0, 'k', 'key', '0.2 m/s');
        set(htv5,'FontSize',16);
    
        [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel(k) vlevel(k)], ...
            'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1 1 1], 'FaceAlpha', 1);
        % c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
        % c2.EdgeColor = [0 0 0];

        m_grid('tickdir','in', ...
       'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
       'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
       'ytick',([38 39]), ... % latitude        
       'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
        
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, string(date));
        % ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
    %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        figname  = [plots_path, 'vort_plots/hslice/seamount_2019_vort'];
        
        figname = strcat(figname, '_vlevel_', string(vlevel(k)));
        figname = strcat(figname, '_nt_', string(indxRange(1) + nt - 1));
        figname = strcat(figname, '.pdf');
        exportgraphics(figure1, figname, 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        % exportgraphics(figure1,strcat(figname, '.eps'))
        close all;
    end
end
% 
%%

%------------------------------------

