function [Var, Att, var_att] = nc2mat(filename)
% nc2mat imports netCDF files (.nc extensions) into MATLAB and stores the 
% variables and attributes contained in the file into 3 structures (.mat)
%
% INPUT:
% filename = 'netCDF_file.nc'
%
% OUTPUT:
% Var = Stucture containing all variables of the .nc file
% Att = Structure containing all attributes of the .nc file
% Var_att = Structure containing the corresponding attributes of each 
% variables
%
% EXAMPLE:
% [Var, Att, var_att] = nc2mat('temperature_Argo_latest.nc')
%
% Created by Earl Duran on 30/7/14
% Reviewed on 1/8/14

nc = ncinfo(filename); % import netCDF file

nc_var_length = length(nc.Variables);
for ii = 1 : nc_var_length
    nc_mat_vers = nc.Variables(1,ii); % select one variable
    var_name_dir = nc_mat_vers.Name; % save its name for the directory
    var_name = var_name_dir; % and again to check the name syntax later
    if isvarname(var_name) == 0
        % changes the name syntax to a valid one if it's not tolerated in
        % MATLAB
        var_name = genvarname(var_name);
    end
    % store the current variable into a structure (using ncread)
    Var.(var_name)=ncread(filename,var_name_dir);
    
    % same method for attributes embedded in the current variable
    var_att_length = length(nc_mat_vers.Attributes);
    if var_att_length > 0
        for jj = 1 : var_att_length
            nc_mat_sub_vers = nc_mat_vers.Attributes(1,jj);
            att_name_dir = nc_mat_sub_vers.Name;
            att_name = att_name_dir;
            if isvarname(att_name) == 0
                att_name = genvarname(att_name);
            end
            % store the current attribute into a structure (using ncreadatt)
            var_att.(var_name).(att_name) = ncreadatt(filename,var_name_dir,att_name_dir);
        end
    end
end

% same method for general attributes
nc_att_length = length(nc.Attributes);
for ii = 1 : nc_att_length
    nc_mat_vers = nc.Attributes(1,ii);
    att_name_dir = nc_mat_vers.Name;
    att_name = att_name_dir;
    if isvarname(att_name) == 0
        att_name = genvarname(att_name);
    end
    Att.(att_name)=ncreadatt(filename,'/',att_name_dir);
end

end