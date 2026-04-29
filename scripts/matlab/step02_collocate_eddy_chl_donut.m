clear all; close all; clc;
%% Paths
repo_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_dir, 'src', 'matlab')));

%% User paths
data_root = '/Users/pearl/Library/CloudStorage/OneDrive-UniversityofTasmania/Work/eac_eddy_env_control_chl_data';

eddy_dir = fullfile(data_root, 'raw', 'aviso_meta');
rs_dir   = fullfile(data_root, 'processed', 'remote_sensing');
out_dir  = fullfile(data_root, 'processed', 'eddy_chl_collocation');

%% Import eddies
eddy_dir = '/Users/renjiongqiu/Work/Documents/eddy_tracking/';
AEdata = nc2mat(strcat(eddy_dir,'anticyclonic_eddy_traj_eac_days_more1_meta3p2_dt_2sat.nc'));
CEdata = nc2mat(strcat(eddy_dir,'cyclonic_eddy_traj_eac_days_more1_meta3p2_dt_2sat.nc'));
% Filtering Virtual Eddies in META dataset
%{
To ensure the integrity of the data, we remove virtual eddies from the Meta dataset, as they may exhibit unusual contours. For this purpose, I've implemented a function named `omitvirtualeddy4mMETA.m`, which should be found in the `utils` directory.
%}
CEvar = omitvirtualeddy4mMETA(CEdata);
AEvar = omitvirtualeddy4mMETA(AEdata);
clear AEdata CEdata
%%
% import chl data
load(fullfile(rs_dir, 'Chl_MODIS_8day_2x4', 'Full_time_Series.mat'));
load(fullfile(rs_dir, 'Chl_MODIS_8day_2x4', 'Chl_fields.mat'));

%% Donut main II
% in this updated version, I stop calculating mld in each contour due to
% the relative coarse resolution of mld dataset. Instead, I calculate raw chl in each contour, 
% for later doing trapping potential analysis. 
dr = 0.2;
nr = 1/dr; % how many dr is in a normalized r=1
r_range = dr:dr:1.5;
donut_struct = struct('CE',struct('chlraw_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'chlraw_dr_velssh',nan(length(CEvar.('time')),1), ...
                                'chlanom_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'chlanom_dr_velssh',nan(length(CEvar.('time')),1), ...
                                'Nobs_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'Nobs_dr_velssh',nan(length(CEvar.('time')),1)), ...
                    'AE',struct('chlraw_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'chlraw_dr_velssh',nan(length(AEvar.('time')),1), ...
                                'chlanom_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'chlanom_dr_velssh',nan(length(AEvar.('time')),1), ...
                                'Nobs_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'Nobs_dr_velssh',nan(length(AEvar.('time')),1)));
tbins_8day = [1:8:365,367];
IE_list = fieldnames(donut_struct);

date_ref = datetime(1950,1,1);
for iyr = 2:length(years)
    yr = years(iyr);
    for ist = 1:length(tbins_8day)-1
        tic
        st = tbins_8day(ist);
        ed = tbins_8day(ist+1);
        % import chl and anomaly field
        iChl = squeeze(Chl(iyr,ist,:,:));
        iChl_loess = squeeze(Chl_loess(iyr,ist,:,:));
        iChl_anom = (iChl - iChl_loess);
        
        % collocate Eddy and chl/anomaly field---------------------------
        for iie = 1:length(IE_list)
            ie = IE_list{iie};
            st_eddy = days( datetime(yr,1,1)+days(st-1) - date_ref );
            ed_eddy = days( datetime(yr,1,1)+days(ed-1) - date_ref );
            eddy = eval([ie 'var']);
            idx_ls = find(eddy.time>=st_eddy&eddy.time<ed_eddy);
            % vel
            lat_vel_ls = eddy.('speed_contour_latitude')(:,idx_ls);
            lon_vel_ls = eddy.('speed_contour_longitude')(:,idx_ls);
            % ssh
            lat_ssh_ls = eddy.('effective_contour_latitude')(:,idx_ls);
            lon_ssh_ls = eddy.('effective_contour_longitude')(:,idx_ls);
            % center
            lon0_ls = eddy.('longitude_max')(idx_ls);
            lat0_ls = eddy.('latitude_max')(idx_ls);

            for iidx = 1:length(idx_ls)
                idx = idx_ls(iidx);
                % vel
                ilat_vel = lat_vel_ls(:,iidx);
                ilon_vel = lon_vel_ls(:,iidx);
                % ssh
                ilat_ssh = lat_ssh_ls(:,iidx);
                ilon_ssh = lon_ssh_ls(:,iidx);
                % center
                ilon0 = lon0_ls(iidx);
                ilat0 = lat0_ls(iidx);
                dlon = (ilon_vel-ilon0)/nr;
                dlat = (ilat_vel-ilat0)/nr;
                
                in_vel0 = false(size(Lat));
                % r_range
                for idr = 1:length(r_range)
                    iilon_vel = ilon0 + dlon*idr;
                    iilat_vel = ilat0 + dlat*idr;
                    in_vel1 = inpolygon(Lon,Lat,iilon_vel,iilat_vel);
                    indonut_vel = logical(double(in_vel1) - double(in_vel0));
                    iNobs = sum(indonut_vel,'all');
                    ichlanom = mean(iChl_anom(indonut_vel),'all','omitnan');
                    ichlraw = mean(iChl(indonut_vel),'all','omitnan');
                    donut_struct.(ie).('chlraw_r_range')(idx,idr) = ichlraw;
                    donut_struct.(ie).('chlanom_r_range')(idx,idr) = ichlanom;
                    donut_struct.(ie).('Nobs_r_range')(idx,idr) = iNobs;
                    in_vel0 = in_vel1;
                end
                % ssh - vel
                in_vel = inpolygon(Lon,Lat,ilon_vel,ilat_vel);
                in_ssh = inpolygon(Lon,Lat,ilon_ssh,ilat_ssh);
                indonut_velssh = logical(double(in_ssh) - double(in_vel));
                iNobs_velssh = sum(indonut_velssh,'all');
                ichlanom_velssh = mean(iChl_anom(indonut_velssh),'all','omitnan');
                ichlraw_velssh = mean(iChl(indonut_velssh),'all','omitnan');
                donut_struct.(ie).('chlraw_dr_velssh')(idx) = ichlraw_velssh;
                donut_struct.(ie).('chlanom_dr_velssh')(idx) = ichlanom_velssh;
                donut_struct.(ie).('Nobs_dr_velssh')(idx) = iNobs_velssh;
            end
        end

        Progres= [iyr,ist];
        fprintf(['iyr: ' num2str(iyr) ' ist: ' num2str(ist) ' yr: ' num2str(yr)])

        % make gif
        
        toc
    end
end
%
fn = ['/Users/renjiongqiu/Work/Documents/aa_remote_sensing/donut_chl8day2x4modis_rawNanom_dr' num2str(dr) '.mat'];
save(fn,'donut_struct','-v7.3');
