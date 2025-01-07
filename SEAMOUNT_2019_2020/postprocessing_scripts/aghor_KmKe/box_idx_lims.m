function [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)

% given which box to average over, define lat-lon idx min-max
if (strcmp(box, 'SW'))
    lat_idx_min = 70;
    lat_idx_max = 90;
    lon_idx_min = 120;
    lon_idx_max = 140;
elseif(strcmp(box, 'S'))
    lat_idx_min = 70;
    lat_idx_max = 90;
    lon_idx_min = 140;
    lon_idx_max = 160;
elseif(strcmp(box, 'SE'))
    lat_idx_min = 70;
    lat_idx_max = 90;
    lon_idx_min = 160;
    lon_idx_max = 180;
elseif(strcmp(box, 'E'))
    lat_idx_min = 90;
    lat_idx_max = 135;
    lon_idx_min = 175;
    lon_idx_max = 190;
elseif(strcmp(box, 'NE'))
    lat_idx_min = 135;
    lat_idx_max = 155;
    lon_idx_min = 160;
    lon_idx_max = 180;
elseif(strcmp(box, 'N'))
    lat_idx_min = 135;
    lat_idx_max = 155;
    lon_idx_min = 140;
    lon_idx_max = 160;
elseif(strcmp(box, 'NW'))
    lat_idx_min = 135;
    lat_idx_max = 155;
    lon_idx_min = 120;
    lon_idx_max = 140;
elseif(strcmp(box, 'W'))
    lat_idx_min = 90;
    lat_idx_max = 135;
    lon_idx_min = 115;
    lon_idx_max = 130;
elseif(strcmp(box, 'big_N'))
    lat_idx_min = 135;
    lat_idx_max = 155;
    lon_idx_min = 120;
    lon_idx_max = 180;
elseif(strcmp(box, 'big_S'))
    lat_idx_min = 70;
    lat_idx_max = 90;
    lon_idx_min = 120;
    lon_idx_max = 180;
end

end