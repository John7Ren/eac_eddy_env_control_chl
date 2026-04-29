function [NewGrid, ilat, ilon, Lat, Lon] = Swath_OSU(WholeGrid, Domain)
%%
%                            Read Me                                      %
% This function will extract a swath specified by the given lat and lon 
%
% This function is specific for the Oregon State Ocean Productivity products
%
% It is based on the description of Lat and Lon described here: http://orca.science.oregonstate.edu/faq01.php
%
% It will check if WholeGrid is the 1080 by 2160 data product, the 2160 by 4320 data product, 
%

%                             Notes                                       %
% This cannot wrap around 180E; if I need to pull a swatch that does, check
% ~/calcs/Swath_new for a workaround using logical indecies


%                             Inputs                                      %
% WholeGrid -> This is the entire spatial grid of data you wish to splice.
%              It can have up to 2 more dimensions, but lat and lon must 
%                    be the last two, respectively
%              You can input:
%              a) the global grid (at 1080 by 2160 resolution)
%              b) the global grid (at 2160 by 4320 resolution)
%              c) Modify section (3.) of this script for a differen input
%                 grid size/location

% Domain    -> This is the desired Domain to splice. 
%              [minlat, maxlat, minlon, maxlon]
%              Enter in Degrees N and Degrees E, respectively

%                             Outputs                                     %
% NewGrid    -> Spliced data, 
%               Returned in as many dimensions as WholeGrid provided 
%                   ! Must be is less or equal to than 4
%
% ilat, ilin -> Indexed output
%               Return vector indecies pulled from WholeGrid
%
% Tlat, Tlon -> Lat and Long grids
%            -> Return length(ilat) x length(ilon) sized matrix of lat lon 
%               grid of NewData

% 1. Set Geographic Bound based on input domain
     minlat = Domain(1);  maxlat = Domain(2);
     minlon = Domain(3);  maxlon = Domain(4);

% 2. Flip Grip 
%    Data in hdf files has lat(1) = +90; We want it to be -90
     WholeGrid = flipud(WholeGrid);

% 2. Check Size of WholeGrid
     gridsize = size(WholeGrid);

% 3. Make Lat and Lon Coordinates 
    if gridsize(1)== 1080
    lon = -180 + (0:2159)*(1/6) + 1/12;
    lat = -90  + (0:1079)*(1/6) + 1/12;
    
    elseif gridsize(1)== 2160
    lon  = -180 + (0:4319)*(1/12) + 1/24;
    lat  = -90 + (0:2159)*(1/12)  + 1/24;
   
    else
        error(['This data input is not the right resoution. It must be one of', ...
              'the two resolution global data sets available on the oregon state ocean productivity site \n'])
    end 


% 3. Find index where lat and lon are in desired domain
     londum = lon(minlon<=lon  & lon<=maxlon);
     ilon   = find(minlon<=lon & lon<=maxlon);

     latdum = lat(minlat<=lat & lat<=maxlat);
     ilat   = find(minlat<=lat & lat<=maxlat);

% 4. Pull from WholeGrid 
     WholeGrid=squeeze(WholeGrid);

     if length(size(WholeGrid))==2
        NewGrid = WholeGrid(ilat,ilon);
     elseif length(size(WholeGrid))==3
        NewGrid = WholeGrid(:,ilat,ilon);
     elseif length(size(WholeGrid))==4
        NewGrid = WholeGrid(:,:,ilat,ilon);
     else
        NewGrid = nan;
     end
 
% 5. Build Lat Lon matricies
     Lon = repmat(lon(ilon), length(ilat),1);
     Lat = repmat(lat(ilat)',1,length(ilon));
%%
end


