clear all; close all; clc;
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
dir = '/Users/renjiongqiu/Work/Documents/aa_remote_sensing/';
% import mld data
% load([dir 'MLD_Hycom_8day_2x4/' 'Full_time_Series']);
% import chl data
load([dir 'Chl_MODIS_8day_2x4/' 'Full_time_Series']);
load([dir 'Chl_MODIS_8day_2x4/' 'Chl_fields']);
%% Donut main
% old version. see the next section for new version.
dr = 0.2;
nr = 1/dr; % how many dr is in a normalized r=1
r_range = dr:dr:1.5;
donut_struct = struct('CE',struct('chl_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'chl_dr_velssh',nan(length(CEvar.('time')),1), ...
                                'mld_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'mld_dr_velssh',nan(length(CEvar.('time')),1), ...
                                'Nobs_r_range',nan(length(CEvar.('time')),length(r_range)), ...
                                'Nobs_dr_velssh',nan(length(CEvar.('time')),1)), ...
                    'AE',struct('chl_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'chl_dr_velssh',nan(length(AEvar.('time')),1), ...
                                'mld_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'mld_dr_velssh',nan(length(AEvar.('time')),1), ...
                                'Nobs_r_range',nan(length(AEvar.('time')),length(r_range)), ...
                                'Nobs_dr_velssh',nan(length(AEvar.('time')),1)));
tbins_8day = [1:8:365,367];
IE_list = fieldnames(donut_struct);

date_ref = datetime(1950,1,1);
for iyr = 2:length(years(1:3))
    yr = years(iyr);
    for ist = 1:length(tbins_8day)-1
        tic
        st = tbins_8day(ist);
        ed = tbins_8day(ist+1);
        % import chl and anomaly field
        iChl = squeeze(Chl(iyr,ist,:,:));
        iChl_loess = squeeze(Chl_loess(iyr,ist,:,:));
        iChl_anom = (iChl - iChl_loess);
        % import mld field
        iMLD = squeeze(MLD(iyr,ist,:,:));
        
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
                    ichl = mean(iChl_anom(indonut_vel),'all','omitnan');
                    imld = mean(iMLD(indonut_vel),'all','omitnan');
                    donut_struct.(ie).('chl_r_range')(idx,idr) = ichl;
                    donut_struct.(ie).('mld_r_range')(idx,idr) = imld;
                    donut_struct.(ie).('Nobs_r_range')(idx,idr) = iNobs;
                    in_vel0 = in_vel1;
                end
                % ssh - vel
                in_vel = inpolygon(Lon,Lat,ilon_vel,ilat_vel);
                in_ssh = inpolygon(Lon,Lat,ilon_ssh,ilat_ssh);
                indonut_velssh = logical(double(in_ssh) - double(in_vel));
                iNobs_velssh = sum(indonut_velssh,'all');
                ichl_velssh = mean(iChl_anom(indonut_velssh),'all','omitnan');
                imld_velssh = mean(iMLD(indonut_velssh),'all','omitnan');
                donut_struct.(ie).('chl_dr_velssh')(idx) = ichl_velssh;
                donut_struct.(ie).('mld_dr_velssh')(idx) = imld_velssh;
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
fn = ['/Users/renjiongqiu/Work/Documents/aa_remote_sensing/donut_chl8day2x4modis_mld8day2x4hycom_dr' num2str(dr) '_test1braket.mat'];
save(fn,'donut_struct','-v7.3');
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

%% experiment on donuting
% create matrix of dots
% x = 145:175;
% y = -20:-1:-50;
x = 162:0.1:165;
y=-44:-0.1:-46;
[X,Y] = meshgrid(x,y);

eddy = CEvar;
i = 1000;
ilat_vel = flip(eddy.('speed_contour_latitude')(:,i));
ilon_vel = flip(eddy.('speed_contour_longitude')(:,i));
ilat_ssh = eddy.('effective_contour_latitude')(:,i);
ilon_ssh = eddy.('effective_contour_longitude')(:,i);
ilon_don = mat2col([ilon_ssh,ilon_vel]);
ilat_don = mat2col([ilat_ssh,ilat_vel]);
u_vel = diff(ilon_vel);
v_vel = diff(ilat_vel);
u_ssh = diff(ilon_ssh);
v_ssh = diff(ilat_ssh);
u_don = diff(ilon_don);
v_don = diff(ilat_don);

in_vel = inpolygon(X,Y,ilon_vel,ilat_vel);
in_don = inpolygon(X,Y,ilon_don,ilat_don);

c1 = 'blue';
c2 = 'red';
close all;
figure('Position',[1200,500,500,400])
hold on;
plot(ilon_vel,ilat_vel,color=c1)
plot(ilon_ssh,ilat_ssh,color=c2)
quiver(ilon_vel(1:end-1),ilat_vel(1:end-1),u_vel,v_vel,color=c1)
quiver(ilon_ssh(1:end-1),ilat_ssh(1:end-1),u_ssh,v_ssh,color=c2)
scatter(X,Y,'x','MarkerEdgeColor','black')
scatter(X(in_vel),Y(in_vel),'bo')
scatter(X(in_don),Y(in_don),'^','MarkerEdgeColor','magenta')

%% experiment on expanding polygon
x = 168:0.1:171;
y=-31:-0.1:-34;
[X,Y] = meshgrid(x,y);

i = 15000;
nr = 1/dr;

ilat_vel = eddy.('speed_contour_latitude')(:,i);
ilon_vel = eddy.('speed_contour_longitude')(:,i);
ilon0 = eddy.('longitude_max')(i);
ilat0 = eddy.('latitude_max')(i);

dlon = (ilon_vel-ilon0)/nr;
dlat = (ilat_vel-ilat0)/nr;
ilon_vel = ilon0 + dlon*idr;
ilat_vel = ilat0 + dlat*idr;

in_vel = inpolygon(X,Y,ilon_vel,ilat_vel);

ilon_vel05 = ilon_vel - (ilon_vel-ilon0)/nr;
ilat_vel05 = ilat_vel - (ilat_vel-ilat0)/nr;
ilon_vel2 = ilon_vel + (ilon_vel-ilon0);
ilat_vel2 = ilat_vel + (ilat_vel-ilat0);
ilat_ssh = eddy.('effective_contour_latitude')(:,i);
ilon_ssh = eddy.('effective_contour_longitude')(:,i);
u_vel = diff(ilon_vel);
v_vel = diff(ilat_vel);
u_ssh = diff(ilon_ssh);
v_ssh = diff(ilat_ssh);
u_vel2 = diff(ilon_vel2);
v_vel2 = diff(ilat_vel2);
u_vel05 = diff(ilon_vel05);
v_vel05 = diff(ilat_vel05);

in_vel = inpolygon(X,Y,ilon_vel,ilat_vel);
in_vel2 = inpolygon(X,Y,ilon_vel2,ilat_vel2);
in_vel05 = inpolygon(X,Y,ilon_vel05,ilat_vel05);
in_ssh = inpolygon(X,Y,ilon_ssh,ilat_ssh);

c1 = 'blue';
c2 = 'red';
c3 = 'green';
c4 = 'magenta';
close all;
figure('Position',[1200,500,500,400])
hold on;
plot(ilon_vel,ilat_vel,color=c1)
plot(ilon_ssh,ilat_ssh,color=c2)
plot(ilon_vel2,ilat_vel2,color=c3)
plot(ilon_vel05,ilat_vel05,color=c4)
quiver(ilon_vel(1:end-1),ilat_vel(1:end-1),u_vel,v_vel,color=c1)
quiver(ilon_ssh(1:end-1),ilat_ssh(1:end-1),u_ssh,v_ssh,color=c2)
quiver(ilon_vel2(1:end-1),ilat_vel2(1:end-1),u_vel2,v_vel2,color=c3)
quiver(ilon_vel05(1:end-1),ilat_vel05(1:end-1),u_vel05,v_vel05,color=c4)
scatter(ilon0,ilat0,'square','filled','MarkerEdgeColor','black','MarkerFaceColor','black')
scatter(X,Y,'x','MarkerEdgeColor','black')
scatter(X(in_vel),Y(in_vel),'bo')
scatter(X(in_ssh),Y(in_ssh),'+','MarkerEdgeColor',c2)
scatter(X(in_vel2),Y(in_vel2),'^','MarkerEdgeColor',c3)
scatter(X(in_vel05),Y(in_vel05),'^','MarkerEdgeColor',c4)

area_vel = polyarea(ilon_vel,ilat_vel);
area_vel2 = polyarea(ilon_vel2,ilat_vel2);