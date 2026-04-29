clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load OSU Remote Sensing Data -- Master  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: Jan 11, 2024 
% By:      Tyler Rohr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Read Me                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script can download and process a collection of the remote senseing products
% from the OSU ocean productivity page. Once you specify which resolution,
% variable and domain you want it will automatically generate a folder for that
% varaiable and save your desired gridded products. You can re-run the script to
% update these products as more data becomes available. You can also
% specify a smaller domain than the entire ocean to process and save.
% Available variables includes all major satellites and NPP models (but not
% CAFE) along with other relevant fields (chl, carbon, growth, mld). See
% Part C of User Input for more info.
%
% Model descriptions can be found here: 
%      https://sites.science.oregonstate.edu/ocean.productivity/vgpm.model.php
%
% Addition variables can be added under part C by adding the appropriate url. Optoions found here: 
%      https://sites.science.oregonstate.edu/ocean.productivity/site.php
%
% Instructions
%     1. Save this script and 'Swath_OSU.m' in the same directory
%     2. Specify what variable and resolution you want to download. More
%        information for each choice is provided by expanding 'Options'
%     3. Specify if you want to save the full gridded time series, a
%        climatology, or both.
%     4. Run Script.
%     5. Re-run script at any point to update your time series. This will
%        look to see if a file already exists then download and append all
%        new data.
%
% List of contents of output save in 'Full_Time_Series' and 'Climatology' output .mat files
%    'dataname' & 'units'      are strings with the respective metadata
%    'Lat' & 'Lon'             are grids of the CENTERED lat and lon of each grid cell
%    'years'                   are the years of data include
%    'time_steps'              are the months of day of year of each time step
%    'RS_data' is the data gridded product product with dimensions of either  
%           -> [years, time_step, Lat(1,:), Lon(:,1)  for the 'Full_Time_Series'
%           -> [time_step, Lat(:,1), Lon(1,:)]        for the 'Climatology'
%    'yticks_m' & 'ylabel_m'   can be used for plotting with monthly ticks.
%                              Note: a table to convert index to day of year is provided below 
%                                    'Save andLabel' at the bottom of the script 
%
% Trouble Shooting
  % 1. Download Urls may be updated after new releases. If they dont work
  %    check the Annual Tar File hyperlinks (right click -> copy link) and
  %    compare to url_var strings assinged in User Inputs Part C.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               User Inputs               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part A. Pick temporal resolution 
  T_Res = 2;
  % Options
   % 1. Monthly Resolution 
   % 2. 8-day Resolution
   for Allocate=1

    if      T_Res == 1
        url_t_res = 'monthly';
       time_steps = 1:12;
    elseif  T_Res == 2 
        url_t_res = '8day';
       time_steps = [1:8:361] + 4; % ~Centered days of averaging 8-day averaging window
                                   %  Remove the '+4' for the starting day
                                
    end

   end

% Part B. Pick spatial resolution 
  S_Res = 2;
  % Options
   % 1. 1/6 degree Resolution  (1080 by 2160)
   % 2. 1/12 degree Resolution (2160 by 4320)
   for Allocate=1

    if     T_Res ==1
       url_s_res = '1x2';
    elseif T_Res == 2 
       url_s_res = '2x4';
    end

   end

 % Part C. Pick variable
  Var = 13;
  % Options
  % 
   % VIIRS Satellite (2012-present) 
     % 1. NPP (mg C/m2/day):        Vertically Generialized Productivity Model (VGPM)
     % 2. NPP (mg C/m2/day): Eppley-Vertically Generialized Productivity Model (eVGPM)
     % 3. NPP (mg C/m2/day):                   Carbon-based Productivity Model (CbPM)
     % 4. Carbon (mg C/m3):                    Carbon-based Productivity Model (CbPM)
     % 5. Growth (1/d):                        Carbon-based Productivity Model (CbPM)
     % 6. Chl    (mg Chl/m3):
     % 7. SST    (C)
   %
   % MODIS Satellite (2002-present) 
     % 8. NPP (mg C/m2/day):        Vertically Generialized Productivity Model (VGPM)
     % 9. NPP (mg C/m2/day): Eppley-Vertically Generialized Productivity Model (eVGPM)
     % 10. NPP (mg C/m2/day):                  Carbon-based Productivity Model (CbPM)
     % 11. Carbon (mg C/m3):                   Carbon-based Productivity Model (CbPM)
     % 12. Growth (1/d):                       Carbon-based Productivity Model (CbPM)
     % 13. Chl    (mg Chl/m3): 
     % 14. SST    (C): 
     % 22. K490    (1/m):
   %
   % SeaWifs Satellite (1997-2009) 
     % 15. NPP (mg C/m2/day):        Vertically Generialized Productivity Model (VGPM)
     % 16. NPP (mg C/m2/day): Eppley-Vertically Generialized Productivity Model (eVGPM)
     % 17. NPP (mg C/m2/day):                   Carbon-based Productivity Model (CbPM)
     % 18. Carbon (mg C/m3):                    Carbon-based Productivity Model (CbPM)
     % 19. Growth (1/d):                        Carbon-based Productivity Model (CbPM)
     % 20. Chl    (mg Chl/m3):
     % XX. SST    (C)                           NOT AVAILABLE HERE because
     %                                          it uses a different grid
   %
   % HYCOM Reanalysis
     % 21. MLD (m) 
   %
   for Allocate=1
       if Var     == 1
       url_var    =  '/vgpm.r2022.v.chl.v.sst/hdf/vgpm.v.';
       long_name  =  'NPP_VGPM_VIIIRS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2012; 
       elseif Var == 2 
       url_var    = '/eppley.r2022.v.chl.v.sst/hdf/eppley.v.';
       long_name  =  'NPP_eVGPM_VIIIRS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2012;
       elseif Var == 3 
       url_var    = '/cbpm2.viirs.r2022/hdf/cbpm.v.';
       long_name  =  'NPP_CbPM_VIIIRS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2012;
       elseif Var == 4 
       url_var    = '/carbon.viirs.r2022/hdf/carbon.v.';
       long_name  =  'Carbon_VIIIRS';
       short_name =  'Carbon';
       units      =  'mgC m$^3$';
       startyear  =  2012;
       elseif Var == 5 
       url_var    = '/growth.viirs.r2022/hdf/growth.v.';
       long_name  =  'Growth_VIIIRS';
       short_name =  'Growth';
       units      =  'd$^{-1}$';
       startyear  =  2012;
       elseif Var == 6 
       url_var    = '/chl.viirs.r2022/hdf/chl.v.';
       long_name  =  'Chl_VIIIRS';
       short_name =  'Chl';
       units      =  'mgChl m$^3$';
       startyear  =  2012;
       elseif Var == 7 
       url_var    = '/sst.viirs.r2022/hdf/sst.v.';
       long_name  =  'SST_VIIIRS';
       short_name =  'SST';
       units      =  '$^\circ$C';
       startyear  =  2012;

       elseif Var == 8
       url_var    =  '/vgpm.r2022.m.chl.m.sst/hdf/vgpm.m.';
       long_name  =  'NPP_VGPM_MODIS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2002;
       elseif Var == 9 
       url_var    = '/eppley.r2022.m.chl.m.sst/hdf/eppley.m.';
       long_name  =  'NPP_eVGPM_MODIS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2002;
       elseif Var == 10 
       url_var    = '/cbpm2.modis.r2022/hdf/cbpm.m.';
       long_name  =  'NPP_CbPM_MODIS';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  2002;
       elseif Var == 11
       url_var    = '/carbon.modis.r2022/hdf/carbon.m.';
       long_name  =  'Carbon_MODIS';
       short_name =  'Carbon';
       units      =  'mgC m$^3$';
       startyear  =  2002;
       elseif Var == 12 
       url_var    = '/growth.modis.r2022/hdf/growth.m.';
       long_name  =  'Growth_MODIS';
       short_name =  'Growth';
       units      =  'd^${-1}$';
       startyear  =  2002;
       elseif Var == 13 
       url_var    = '/chl.modis.r2022/hdf/chl.m.';
       long_name  =  'Chl_MODIS';
       short_name =  'Chl';
       units      =  'mgChl m$^3$';
       startyear  =  2002;
       elseif Var == 14 
       url_var    = '/sst.modis.r2022/hdf/sst.m.';
       long_name  =  'SST_MODIS';
       short_name =  'SST';
       units      =  '$^\circ$C';
       startyear  =  2002;
       elseif Var == 22 
       url_var    = '/k490.modis.r2022/hdf/k490.m.';
       long_name  =  'K490_MODIS';
       short_name =  'K490';
       units      =  '$m^{-1}$';
       startyear  =  2002;
       elseif Var == 23 
       url_var    = '/par.modis.r2022/hdf/par.m.';
       long_name  =  'PAR_MODIS';
       short_name =  'PAR';
       units      =  '$m^{-2}$ d$^{-1}$';
       startyear  =  2002;

       elseif Var == 15
       url_var    =  '/vgpm.r2022.s.chl.a.sst/hdf/vgpm.s.';
       long_name  =  'NPP_VGPM_Seawifs';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  1997;
       elseif Var == 16 
       url_var    = '/eppley.r2022.s.chl.a.sst/hdf/eppley.s.';
       long_name  =  'NPP_eVGPM_Seawifs';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  1997;
       elseif Var == 17 
       url_var    = '/cbpm2.seawifs.r2022/hdf/cbpm.s.';
       long_name  =  'NPP_CbPM_Seawifs';
       short_name =  'NPP';
       units      =  'mgC m$^2$ d$^{-1}$';
       startyear  =  1997;
       elseif Var == 18
       url_var    = '/carbon.seawifs.r2022/hdf/carbon.s.';
       long_name  =  'Carbon_Seawifs';
       short_name =  'Carbon';
       units      =  'mgC m$^3$';
       startyear  =  1997;
       elseif Var == 19
       url_var    = '/growth.seawifs.r2022/hdf/growth.s.';
       long_name  =  'Growth_Seawifs';
       short_name =  'Growth';
       units      =  'd$^{-1}$';
       startyear  =  1997;
       elseif Var == 20
       url_var    = '/chl.seawifs.r2022/hdf/chl.s.';
       long_name  =  'Chl_Seawifs';
       short_name =  'Chl';
       units      =  'mgChl m$^3$';
       startyear  =  1997;
       % elseif Var == 21
       % url_var    = '/sst.avhrr/hdf/sst.a.';
       % long_name  =  'SST_Seawifs';
       % short_name =  'SST';
       % units      =  '$^\circ$C';
       % startyear  =  1997;

       elseif Var == 21
       url_var    = '/mld125.hycom/hdf/mld.hycom_125.';
       long_name  =  'MLD_Hycom';
       short_name =  'MLD';
       units      =  '$m$';
       startyear  =  1997;


    end

   end

  % Part D. Pick domain
   Domain_Range = 2; 
   % Options
    % 1. Global: Downloads full global files
    % 2. Custom: Splices out specified Lat/Lon Box for smaller file
    %            storage. User must specify under Allocate tab below
    for Allocate=1

    if Domain_Range == 1
       Min_lat = -90;  Max_lat = 90;   % Degrees N (so S is negative)
       Min_lon = -180;  Max_lon = 180; % Degrees E (so W is negative)
       Domain = [Min_lat,Max_lat,Min_lon,Max_lon];
    elseif Domain_Range == 2 
       % Update:  As Desired.
       %   NOTE: Domain cannot cross 180E without modifying the Swath_OSU script. 
       Min_lat = -50;  Max_lat = -20;   % Degrees N (so S is negative)
       Min_lon = 145;  Max_lon = 175; % Degrees E (so W is negative)
       Domain = [Min_lat,Max_lat,Min_lon,Max_lon];
    end

   end
  
 % Part E. Pick time range
   Time_Range = 2; 
   % Options
    % 1. All Available: Downloads from Satelittle Launch data (see variable
    %                   options (Part C) until present day.
    % 2. Custom:        Downloads specified year range. User must specify
    %                   under Allocate tab below.  !NOTE! Using this option
    %                   will overwrite any previous file if a new root diretory is not prescribed 
  for Allocate=1

    if     Time_Range == 1
          % startyear          is assigened automatically based on the
                             % variable/satelitte chosen in Part C
            lastyear = 2030; % The script will try and download all files up to this year 
                             % but will not do anything for years with no
                             % files. Will need to update after 2030.     
    elseif Time_Range == 2 
           startyear = 2002;
           lastyear  = 2023;
 
    end

  end

 % Part E. Pick root directory to save output & directory where Swath_OSU.m is saved
   root_directory = '/Users/pearl/Library/CloudStorage/OneDrive-UniversityofTasmania/Work/data';                 % TRs path on casper
   
 % Part F. Pick output to save
   Save_Full_Time_Series = true;
   Save_Full_Climatology = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Load Data              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   
% 1. Create directory and file names 
     data_dir  = strcat(root_directory,'/',long_name,'_',url_t_res,'_',url_s_res );
     dataname = strcat(long_name,'_',url_t_res,'_',url_s_res );
     addpath(pwd)

% 2. Check for existing directory 
     if exist(data_dir, 'dir') == 7 & Time_Range == 1
        fprintf([ '\n', dataname, ' already exists and the full time domain has been requested.', ...
                  'Appending any new data and updating the climatology. \n'] );
        update_flag = true;
        cd(data_dir)
        load Full_Time_Series
        % Start with last year included (and eventually replace) to fill in any new days/months added 
          startyear = years(end); years_old = years; clearvars years;
        
        % Loop through years starting at one year before the last year in
        % the existing file
        for y =startyear:lastyear 
            % Write hyperlinked url
            url = strcat('http://orca.science.oregonstate.edu/data/', url_s_res,'/', url_t_res,  url_var, num2str(y),'.tar');
            % download to directory
            try
            % download to directory if link exists
            websave([data_dir,'/',short_name,num2str(y)], url);
            fprintf(['\n Donwloaded .tar file for year: ' num2str(y)])
            catch
            % delete html error file if link does not exists
            delete *html
            end
        end


     else
        fprintf([ '\n', dataname, ' has not yet been downloaded. Downloading the full time series. \n'] );
        update_flag = false;
        mkdir(data_dir); cd(data_dir)

        % Loop through years starting at first year of satellite coverage until files dont exist anymore
        for y =startyear:lastyear 
            % Write hyperlinked url
            url = strcat('http://orca.science.oregonstate.edu/data/', url_s_res,'/', url_t_res,  url_var, num2str(y),'.tar');
            fprintf(['looping ...' num2str(y) '...'])
            % download to directory
            try
            % download to directory if link exists
            websave([data_dir,'/',short_name,num2str(y)], url);
            fprintf(['\n Susscuflly Donwloaded Year: ' num2str(y)])
            catch
            % delete html error file if link does not exists
            delete *html
            end
        end

     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Process Data              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 2. Compile list of all annual (.tar) files
     dir_list_tar=dir('*.tar'); 
     dir_names_tar={dir_list_tar(:).name};

% 3. Loop through each year and build matlab structure

     for i = 1:length(dir_list_tar)
         % Extract year number from file name and print progress
           dumSTR = dir_names_tar{i}; y_str = dumSTR(end-7:end-4); years(i)  =  str2num(y_str);  
           fprintf([' \n \n .... Processing Year: ' y_str ' .... \n'])    
    
          % Untar and unzip file
            try
            untar(dir_names_tar{i});  % releases zipped daily hdf files
            catch
            fprintf(['Problem with .tar file, skipped year:', y_str, '\n']);
            end
          % Unzips all hdf files
            gunzip('*.gz')             
            
          % List all hdf files names in folder
            dir_list_hdf=dir('*.hdf'); 
            dir_names_hdf={dir_list_hdf(:).name};
    
          % Preallocate matlab structure 
                   
                     if i==1  % (only run on first loop)
            
                       % Time Steps per year
                         ts_length = length(time_steps);  
            
                       % Years 
                         years_length = length(dir_names_tar); 
                
                       % Spatial Domain
                       % Load dummy data to get coordinate indecies
                          fname    = dir_names_hdf{1};        % Set filename to first day of data
                          hinfo    = hdfinfo(fname);          % Pull metadata from first day of data
                          varname  = hinfo.SDS.Name;          % Pull variable name
                          data_dum = hdfread(fname, varname); % Pull data
                          
                        % Splice Domain
                        [data, ilat, ilon, Lat, Lon]  =  Swath_OSU(data_dum, Domain);  
                        
                        % Preallocate Matrix  
                        DATA = nan(years_length,ts_length, length(ilat), length(ilon));
                        fprintf('\n Matrix Preallocatted. \n')
                     end
 
          % Loop through daily hdf files and add to SOTS_DATA matrix 

                    for ts = 1:length(dir_names_hdf)
                        fprintf(['\n Processing timestep = ' num2str(ts) '\n'])
                    
                        % Find time step index based on explict name of .hdf file
                            file_name = dir_names_hdf{ts};
                            % Monthly
                            if T_Res == 1
                            Start_day=[1 32	61 92 122 153 183 214 245 275 305 335]; 
                            doy_index = (str2num(file_name(end-6:end-4))); 
                            [dum,ts_index] = (min(abs(doy_index - Start_day)));
                            elseif T_Res == 2
                            ts_index =(str2num(file_name(end-6:end-4)) + 7)/8; 
                            end
            
                        % Pull variable data from .hdf file
                            data_dum = hdfread(file_name,varname);
                         
                        % Pull desired swath out of global field
                           if     Domain_Range == 1
                                  % No need to splice global domain
                                  data=data_dum;
                           elseif Domain_Range == 2
                                  % Splice out desired Section
                                  [data, ilat, ilon, TLat, TLon]  =  Swath_OSU(data_dum, Domain);
                           end
                      
                        % Put in DATA Matrix
                            DATA(i,ts_index, :,:) =  data;
                    end
        
        
        
           % Delete .hdf and .tar files
              delete *.hdf*
     end

% 4. Delete .tar files
      delete *.tar 

% 5. Remove Land and replace with nan
      fprintf('\n Replacing land with nan \n')
      Land_ind = find(DATA==-9999 | DATA == -833.2500); 
      DATA(Land_ind) = nan; 
%% check if variables exist:
% update_flag, dataname
% 6. Add to old matrix if it exists
      if update_flag
          fprintf('\n Appending new data onto prior matrix \n')
          eval(['Data_old = ', short_name ,';'])
          DATA = cat(1,Data_old(1:end-1,:,:,:),DATA);
          years = cat(2,years_old(1:end-1),years);
      end

% 7. Create climatology
     DATA_Clim=squeeze(mean(DATA,'omitnan'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Save and Label              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make Timing lables. Expand for Table on day of year conversions
%    Month     (doy)        % indecies    % Corresponding centered days   (doys)
%    January   (1:30)   :       1-4           Jan. 5th - Jan. 29          (5-29)
%    Februrary (32:59)  :       5-7           Feb. 6th - Feb. 22          (37-53)    
%    March     (60:90)  :       8-11          Mar. 2th - Mar. 26          (61-85)
%    April     (91:120) :       12-15         Apr. 3rd - Apr. 27          (93-117)
%    May       (121:151):       16-19         May. 5th - May. 29          (125-149)
%    June      (152:181):       20-23         Jun. 6th - Jun. 30          (157-181)
%    July      (182:212):       24-26         Jul. 8th - Jul. 24          (189-205)
%    August    (213:243):       27-30         Aug. 1th - Aug. 25          (213-237) 
%    September (244:273):       31-34         Sep. 2th - Sep. 26          (245-269) 
%    October   (274:304):       35-38         Oct. 4th - Oct. 28          (277-301)
%    November  (305:334):       39-42         Nov. 5th - Nov. 29          (309-333)
%    December  (335:365):       43-46         Dec. 7th - Dec. 28          (341-365)

     if T_Res == 1
      yticks_m = 1:12;
      ylabel_m = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
  
     elseif T_Res == 2
      yticks_m = [    1      5        8      12      16      20     24      27      31     35       39     43];
      ylabel_m = {'Jan 5' 'Feb 6' 'Mar 2' 'Apr 3' 'May 5' 'Jun 6' 'Jul 8' 'Aug 1' 'Sep 2' 'Oct 4' 'Nov 5' 'Dec 7'};
     end 

% Save

      fprintf('\n Saving: \n')
      
      if Save_Full_Time_Series
      fprintf('\n .... Full Time Series .... \n')
      eval([short_name ,'= DATA;'])
      save('Full_Time_Series' , 'dataname', 'units', short_name, 'Lat', 'Lon', ... 
                                'time_steps', 'years', 'yticks_m','ylabel_m', 'Domain', '-v7.3')
      end

      if Save_Full_Climatology
      fprintf('\n .... Climatology .... \n')
      eval([short_name ,'= DATA_Clim;'])
      save('Climatology',      'dataname', 'units', short_name, 'Lat', 'Lon', ...
                               'time_steps', 'years', 'ylabel_m', 'Domain', '-v7.3')
      end

      fprintf(['\n' dataname ' All done! \n'])

%% -------------------- calculate the loess filter
clear all; close all; clc;
dir = '/Users/renjiongqiu/Work/Documents/aa_remote_sensing/Chl_MODIS_8day_2x4/';
load([dir 'Full_time_Series']);

gLat = imresize(Lat,0.5);
gLon = imresize(Lon,0.5);
Chl_squeeze = nan(length(years)*length(time_steps),360,360);
Chl_loess = zeros(size(Chl));
Chl_32tmean = nan(size(Chl));
Chl_32tmean_loess = zeros(size(Chl));
% 32-d time window
for iyr = 1:length(years)
    for idt = 1:length(time_steps)
        ii = length(time_steps)*(iyr-1) + idt;
        Chl_squeeze(ii,:,:) = Chl(iyr,idt,:,:);
        % disp(ii);
    end
end
Chl_32tmean_squeeze = movmean(Chl_squeeze,5,1,'omitnan');
for iyr = 1:length(years)
    for idt = 1:length(time_steps)
        ii = length(time_steps)*(iyr-1) + idt;
        Chl_32tmean(iyr,idt,:,:) = Chl_squeeze(ii,:,:);
        % disp(ii);
    end
end

for iyr = 1:length(years)
    yr = years(iyr);
    for idt = 1:length(time_steps)
        tic
        dt = time_steps(idt);
        % calculate chl field
        iChl = squeeze(Chl(iyr,idt,:,:));
        giChl = imresize(iChl,0.5);
        [lp,~]=smooth2d_loess(giChl,gLon(1,:),gLat(:,1),6,6,gLon(1,:),gLat(:,1)); %smooth using Loess filter
        lp_data = imresize(lp,1/0.5);
        Chl_loess(iyr,idt,:,:) = lp_data;
        % calculate 32tmeaned chl field
        iChl = squeeze(Chl_32tmean(iyr,idt,:,:));
        giChl = imresize(iChl,0.5);
        [lp,~]=smooth2d_loess(giChl,gLon(1,:),gLat(:,1),6,6,gLon(1,:),gLat(:,1)); %smooth using Loess filter
        lp_data = imresize(lp,1/0.5);
        Chl_32tmean_loess(iyr,idt,:,:) = lp_data;

        Progres= [iyr,idt];
        fprintf(['iyr: ' num2str(iyr) ' idt: ' num2str(idt) ' yr: ' num2str(yr)])
        toc
    end
end
% save('Chl_Anom_TS' , 'dataname', 'units', short_name, 'Lat', 'Lon', ... 
%                                 'time_steps', 'years', 'yticks_m','ylabel_m', 'Domain', '-v7.3')
%
fn = '/Users/renjiongqiu/Work/Documents/aa_remote_sensing/Chl_MODIS_8day_2x4/Chl_fields.mat';
save(fn,'Chl_loess','Chl_32tmean_loess','-v7.3');
