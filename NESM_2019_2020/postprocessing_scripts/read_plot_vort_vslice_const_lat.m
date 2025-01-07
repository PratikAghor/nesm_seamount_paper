%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1

%%Loop to calculate averages if need (modify as needed)
indxRange = 1691:3382;
% indxRange = 2144:2183; % What time indices do you need?
nt0=indxRange(1);
% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

tskip = 1;
vlevel_skip = 50;
%-------------------------------------------------------
%-------------------------------------------------------
%%
vort_t_filename = strcat('nesm_2019_2020_vort_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% vort_t_filename = 'nesm_2019_vort_t_1000_2902.nc';
% nt0 = 1000; % 0 if 0_999, 1000 if 1000_2902
% nt0=indxRange(1);
t_arr = ncread(vort_t_filename, 't_arr');
lat_arr = ncread(vort_t_filename, 'lat');
vort_t_vslice = ncread(vort_t_filename, 'vort_t_vslice');
%%
%-------------------------------------------------------

% mkdir vort_plots;
%% plotting
for nt = 1300:tskip:Nt
    sprintf('for loop indx = %d', nt)
    vort_t_nt = indxRange(1) - nt0 + nt;
    sprintf('vort_t arr indx = %d', vort_t_nt)
    sprintf('fig indx = %d', vort_t_nt + nt0 - 1)
    for k = 2:2 % length(lat_arr)
        f0 = 2*Omega*sin(lat_arr(k)); % Coriolis freq.
        
        lat_idx = find(abs(lat_rho_vec-lat_arr(k))<1e-3); % find idx of lat in lat_rho_vec
        z_depth_vec = squeeze(zr(:, lat_idx, :)); % get depth vals at that lat
        Z = ((squeeze(zr(:,lat_idx,:))));
        X = (repmat(lon_rho_vec, NumLayers, 1));
        %%% plotting
        figure1 = figure(nt);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        vort_vslice = squeeze(vort_t_vslice(vort_t_nt, k, :, :));
        zMin = -0.5; % min(min(vort_vslice)); 
        zMax = -zMin; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        pcolor(ax1, X, Z, transpose(vort_vslice)); 
        shading interp;
        set(ax1,'Color', [1 1 1])
        % set(ax1,'YDir','reverse')
        % colormap(ax1,b2r(zMin,zMax));  
        % colormap(ax1, whitejet); clim([zMin zMax]);
        % colormap(ax1, "jet"); 
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        % clim([zMin zMax]);
        colorbar; 
       
        % title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, string(date));

        ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center

        xticks([-64.99 -64 -63 -62 -61])  % longitude   
        xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
        yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
        ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        ylim([-5000 0])
    %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        figname  = [plots_path, 'vort_plots/vslice/const_lat/nesm_2019_2020_vort_vslice'];
        
        figname = strcat(figname, '_lat_', string(lat_arr(k)));
        figname = strcat(figname, '_nt_', string(indxRange(1) + nt - 1));
        % vort_contour = strcat(vort_contour, '.pdf');
        exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        % set(figure1, 'PaperPositionMode', 'auto')
        % print(figure1,strcat(figname, '.png'),'-dpng','-r300');
        % exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
        % export_fig(ax1, figname, '-eps','-transparent', '-r300');
        close all;
    end
end
% 
%%

%------------------------------------

