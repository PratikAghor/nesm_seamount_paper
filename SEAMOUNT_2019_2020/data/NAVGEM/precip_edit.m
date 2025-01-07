% edit precip file to match frc_time and blk_time 
% before matlab, copy precip_file to edited_precip_file.
% now we will only edit MT, Date and precip attributes in the new file. 
clear all; clc;

% workdir='~/scratch/work/SEAMOUNT_2019_2020/data/NAVGEM/'
precip_file = 'navgem1.4_0.281c-std_2020_03hr_precip.nc'
edited_precip_file = 'navgem1.4_0.281c-std_2020_03hr_precip_edited.nc'

Date_uv = ncread('navgem1.4_0.281c-sec_2020_03hr_uv-10m.nc', 'Date'); % Dates we want
MT_uv = ncread('navgem1.4_0.281c-sec_2020_03hr_uv-10m.nc', 'MT'); % MT we want
Nt_uv=size(MT_uv)

MT_precip=ncread(precip_file, 'MT'); % MT we have
size(MT_precip)

Lat = ncread(precip_file, 'Latitude');
Lon = ncread(precip_file, 'Longitude');
% read from precip 2020 file, requires multiple processors
% precip_old=ncread(precip_file, 'precip'); % does not work

ncid_old = netcdf.open(precip_file,'NC_NOWRITE');
varid_old = netcdf.inqVarID(ncid_old, 'precip');

precip_old = netcdf.getVar(ncid_old, varid_old);

netcdf.close(ncid_old);

[Nlon, Nlat, Nt_precip] = size(precip_old)

precip_old(7, 9, 100)

precip = zeros(Nlon, Nlat, length(MT_uv));

ncid = netcdf.create(edited_precip_file, 'NETCDF4');
netcdf.close(ncid);

ncid = netcdf.open(edited_precip_file,'WRITE');

% define variables and attributes
netcdf.reDef(ncid);

	% 1. MT
	MT_dimid = netcdf.defDim(ncid,'MT', length(MT_uv));
	MT_varid = netcdf.defVar(ncid, 'MT', 'NC_DOUBLE', [MT_dimid]);
	netcdf.putAtt(ncid, MT_varid, 'long_name', 'time');
	netcdf.putAtt(ncid, MT_varid, 'units', 'days since 1900-12-31 00:00:00');
	netcdf.putAtt(ncid, MT_varid, 'calendar', 'standard');
	netcdf.putAtt(ncid, MT_varid, 'axis', 'T');
	netcdf.putAtt(ncid, MT_varid, 'next_MT', 43600.625);

	% 2. Date

	Date_dimid = netcdf.defDim(ncid,'Date', length(Date_uv));
	Date_varid = netcdf.defVar(ncid, 'Date', 'NC_DOUBLE', [Date_dimid]);
	netcdf.putAtt(ncid, Date_varid, 'long_name', 'Date');
	netcdf.putAtt(ncid, Date_varid, 'units', 'day as %Y%m%d.%f');
	netcdf.putAtt(ncid, Date_varid, 'C_format', '%13.4f');
	netcdf.putAtt(ncid, Date_varid, 'FORTRAN_format', '(f13.4)');
	netcdf.putAtt(ncid, Date_varid, 'next_Date', 20200515.625);

	% 3. Lat
	Lat_dimid = netcdf.defDim(ncid,'Latitude', length(Lat));
	Lat_varid = netcdf.defVar(ncid, 'Latitude', 'NC_DOUBLE', [Lat_dimid]);
        netcdf.putAtt(ncid, Lat_varid, 'standard_name', 'latitude');
        netcdf.putAtt(ncid, Lat_varid, 'units', 'degrees_north');
        netcdf.putAtt(ncid, Lat_varid, 'axis', 'Y');
	
	% 4. Lon
	Lon_dimid = netcdf.defDim(ncid,'Longitude', length(Lon));
        Lon_varid = netcdf.defVar(ncid, 'Longitude', 'NC_DOUBLE', [Lon_dimid]);
        netcdf.putAtt(ncid, Lon_varid, 'standard_name', 'longitude');
        netcdf.putAtt(ncid, Lon_varid, 'units', 'degrees_east');
	netcdf.putAtt(ncid, Lon_varid, 'point_spacing', 'even');
	netcdf.putAtt(ncid, Lon_varid, 'modulo', '360 degrees');
        netcdf.putAtt(ncid, Lon_varid, 'axis', 'X');

	% 5. precip
	precip_varid = netcdf.defVar(ncid, 'precip', 'NC_SHORT', [Lon_dimid Lat_dimid MT_dimid]);
	netcdf.putAtt(ncid, precip_varid, 'coordinates', 'Date');
	netcdf.putAtt(ncid, precip_varid, 'time_average_bounds', '[0 0]');
	netcdf.putAtt(ncid, precip_varid, 'long_name', ' precipitation'); % copied from precip_file attributes to the extra space!
	netcdf.putAtt(ncid, precip_varid, 'standard_name', 'lwe_precipitation_rate');
	netcdf.putAtt(ncid, precip_varid, 'units', 'm/s');
	netcdf.putAtt(ncid, precip_varid, 'scale_factor', 7.6294e-06);
	netcdf.putAtt(ncid, precip_varid, 'add_offset', 0.25);

% close define mode
netcdf.endDef(ncid);

% write data
netcdf.putVar(ncid, MT_varid, MT_uv);
netcdf.putVar(ncid, Date_varid, Date_uv);
netcdf.putVar(ncid, Lat_varid, Lat);
netcdf.putVar(ncid, Lon_varid, Lon);

% now interpolate/extrapolate precip vals
for i=1:Nlon
	for j=1:Nlat
		for k=1:Nt_precip
			precip(i, j, k) = precip_old(i, j, k);
		end
	end
end

disp('precip_old, precip')
precip_old(7, 9, 100)
precip(7,9, 100)
% extrapolate
for i=1:Nlon
	for j=1:Nlat
		for k=(Nt_precip + 1):length(MT_uv)
			precip(i, j, k) = 2*precip(i, j, k-1) - precip(i, j, k-2); 
		end
	end
end

netcdf.putVar(ncid, precip_varid, precip);
netcdf.close(ncid);

clear all
