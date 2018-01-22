% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%Finds GHCN data for coordinates, dates, and metVar specified and writes to
%ASCII file.

crd = [-72, -64, -26, -17];
dates = [1970, 2012];
region = 'High_Andes';


%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules


%%INNERWORKINGS (DO NOT MODIFY):
rootPathData = fullfile(pwd, 'GHCN_data');

for ii = 2 : -1 : 1
    if ii == 1
        dataCon = 'adj';
    elseif ii == 2
        dataCon = 'non';
    end

    for jj = 1 : 2
        if jj == 1
            metVar = 'pre';
        elseif jj == 2
            metVar = 'tmn';
        end
        
        useData = GHCN_find(crd, dates, metVar, dataCon);

        fullPathData = fullfile(rootPathData, [region, '_GHCN_', metVar, '_', dataCon, '.asc']);
        write_GHCN(useData, fullPathData);
    end
end

