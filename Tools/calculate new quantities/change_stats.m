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

clear all

%%SELECT MODE:
inputs = 'abs_ensemble'; 
    %'abs_ensemble': will be promped for (1) ensemble projection grid, 
        %(2) standard deviation of ensemble, and (3) hindcast grid.
    %'change_ensemble': will be promped for (1) ensemble change grid and
        %(2) standard deviation of ensemble
    %'abs_single': will be promped for (1) projection grid from one model 
        %and (2) hindcast grid.
        
startPath = 'E:\Figure Repository\MHS_GCM_2014 article figures\aggregated change table\data';


%GUI to get files:
if ~exist('startPath','var') || ~exist(startPath,'dir')
    startPath = pwd;
end

if regexpbl(inputs,'abs')
    [fileData, foldData, ~] = uigetfile({'*.asc'}, ...
        'Select projected quantity', startPath);
    disp([fullfile(foldData,fileData) ' was selected as projected quantity.']);
elseif regexpbl(inputs,'change')
        [fileData, foldData, ~] = uigetfile({'*.asc'}, ...
        'Select change in projected quantity', startPath);
    disp([fullfile(foldData,fileData) ' was selected as change in projected quantity.'])
end

if regexpbl(inputs,'ensemble')
    [fileStd, foldStd, ~] = uigetfile({'*.asc'}, ...
        ['Select ensemble' char(39) 's standard deviation of projected quantity'], foldData);
    disp([fullfile(foldStd,fileStd) ' was selected as the ensemble standard deviation.'])
end

if regexpbl(inputs,'abs')
   indSrch = regexpi(foldData,filesep);
    pathSrch = foldData(1:indSrch(end-2));
    [fileHis, foldHis, ~] = uigetfile({'*.asc'}, ...
        'Select hindcast data asscoiated with projected quantity', pathSrch);  
    disp([fullfile(foldHis,fileHis) ' was selected as the hindcast data.'])
end

%%LOAD DATA:
[data, hdrData, metaData] = read_ESRI(fullfile(foldData,fileData));

if regexpbl(inputs,'abs')
    [hind, hdrHis, metaHis] = read_ESRI(fullfile(foldHis,fileHis));
    
    if ~isequal(hdrData,hdrHis)
        error('agg_change:misAlign','The projected grid and associated hindcast grids do not align.');
    end
    
    %Calculate Change:
    change = data - hind;
else
    change = data;
end

if regexpbl(inputs,'ensemble')
    [sd, hdrSd, metaSd] = read_ESRI(fullfile(foldStd,fileStd));
    
    %Error if grids don't align:
    if ~isequal(hdrData,hdrSd)
        error('agg_change:misAlign','The projected and associated SD grids do not align.');
    end
    
end


%Calculate area of region
[lat, lon] = ESRI_hdr2geo(hdrData,metaData);
area = area_geodata(lon, lat, 'c');

%%DISPLAY RESULTS:
if regexpbl(inputs,'ensemble')
    %Calculate one SD away from mean change:
    upper = sum2d((change+sd).*area)/10^3;
    totChange = sum2d(change.*area)/10^3;
    lower = sum2d((change-sd).*area)/10^3;

    aggSd = mean(upper-totChange, totChange - lower);
    %Display results:
    disp(['Mean = ' num2str(totChange,'%.2E') 'm^3 (' num2str(lower,'%.2E') ' to ' num2str(upper,'%.2E') '; SD = \pm ' num2str(aggSd,'%.2E') ')']);
    disp(['Average loss per area = ' num2str(totChange/sum2d(area),'%.2E') 'm^3 (' num2str(lower/sum2d(area),'%.2E') ' to ' num2str(upper/sum2d(area),'%.2E') ')']);
else
    %Display results:
    disp(['Mean = ' num2str(totChange,'%.2E') 'm^3']);
    disp(['Average loss per area = ' num2str(totChange/sum2d(area),'%.2E') 'm^3']);
end
