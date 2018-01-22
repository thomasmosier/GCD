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

%Script finds info from ASCII files in selected folder and tabulated info
%based upon defined queries:


%If selected (1), all snow time-series grids are written as geospatial ASCII files.
wrtTs = 1;
fileTyp = 'ascii'; %Choose 'netCdf' or 'ascii'

%Only need these if NetCDF chosen:
yrs = [2020,2100]; %Only used if netCDF chosen.
crd = [64, 86, 30.5, 39]; %Only used if netCDF chosen.

%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end-2); %Assumes scripts are in distribution hierarchy.

addpath(genpath(pathScript(1:indParent-1))); %Add modules

if ~isempty(regexpi(fileTyp,'asc'))
    %Get folder with ASCII temp data:
    strTempGet = ['Select the directory containing ESRI ASCII monthly '...
        'temperature data'];
    foldTmp = uigetdir(pwd,strTempGet);
    disp([foldTmp ' has been selected as the folder of monthly temperature files.']);

    %Get folder with ASCII precip data:
    strPreGet = ['Select the directory containing ESRI ASCII monthly '...
        'precipitation data'];
    indTmpSep = regexpi(foldTmp,filesep);
    dirStartPre = foldTmp(1:indTmpSep(end-1)-1);
    foldPre = uigetdir(dirStartPre, strPreGet);
    disp([foldPre ' has been selected as the folder of monthly Precipitation files.']);
    
elseif ~isempty(regexpi(fileTyp,'cdf'))
    strTempGet = ['Select the NetCDF file containing ESRI ASCII monthly '...
    'temperature data'];
    [fileTmp, foldTmp] = uigetfile({'*.nc';'*.grb'},strTempGet,pwd);
    fileTmp = fullfile(foldTmp,fileTmp);
    disp([fileTmp ' has been selected as the monthly temperature source.']);

    indTmpSep = regexpi(foldTmp,filesep);
    dirStartPre = foldTmp(1:indTmpSep(end-1)-1);
    
    strPreGet = ['Select the NetCDF file containing ESRI ASCII monthly '...
    'precipitation data'];

    [filePre, foldPre] = uigetfile({'*.nc';'*.grb'},strPreGet,dirStartPre);
    filePre = fullfile(foldPre,filePre);
    disp([filePre ' has been selected as the monthly precipitation source.']);
else
    error('snowfall_calc:fileType',[fileTyp ' is an unknown file type, which has not been coded for.'])
end

if regexpbl(fileTyp,'asc')
    pathSnow = snowfall(foldPre,foldTmp,fileTyp,wrtTs);
elseif regexpbl(fileTyp,'cdf')
    pathSnow = snowfall(foldPre,foldTmp,fileTyp,wrtTs,crd,yrs);
end