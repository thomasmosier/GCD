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

function [dataTsShep, corrMat, corrMeta, stnData] = shepard_corr_w_GHCN(sTsScaled,sMeta,param,varargin)


p = param(1);
r = param(2);

%Geographic extent used to find GHCN stations can be set larger than
%downscaled grid extent.
%     GHCNcrd(1) = crd(1) - hdrTsRaw(5)/fracLap;
%     GHCNcrd(2) = crd(2) + hdrTsRaw(5)/fracLap;
%     GHCNcrd(3) = crd(3) - hdrTsRaw(5)/fracLap;
%     GHCNcrd(4) = crd(4) + hdrTsRaw(5)/fracLap;

if length(size(sTsScaled.data)) == 2
    currTs = single(sTsScaled.data);
elseif length(size(sTsScaled.data)) == 3
    [indTs, ~] = geodata_time_ind(sTsScaled, sMeta.currTime);
    currTs = single(squeeze(sTsScaled.data(indTs,:,:)));
else
    error('delta_calc:unkownDim',['The dimenstions of the current '...
        'array are ' num2str(length(size(sLrTs.data)))] );
end

dataTsShep = nan(size(currTs),'single');

if length(varargin) >= 2 && ~isempty(varargin{1})
    stnData = varargin{1};
    strGhcnTyp = varargin{2};
elseif length(varargin) == 0 || isempty(varargin{1}) %#ok<ISMT>
    %Find corresponding GHCN data:
    stnData = GHCN_find(sMeta.crd, sMeta.yrsOut, sMeta.currVar, 'adj');
        %stnData is a cell array where: 
        %{i,1} is the GHCN station ID,
        %{i,2} is the station latitude, 
        %{i,3} is the station longitude,
        %{i,4} is the station elevation, and 
        %{i,5} is an m by 13 cell-array with yearly data (first column is year).
        
    strGhcnTyp = 'adj';    
    if isempty(stnData(:,1))    %If no adjusted GHCN stn records, try finding non-adjusted records.
        disp(['No adjusted GHCN stations were found for ' sMeta.region '.']);

        stnData = GHCN_find(sMeta.crd, sMeta.yrsOut, sMeta.currVar, 'org');
            warning('grid_correct:nonAdjust', ['Non-adjusted '...
                'GHCN station records are being used in the '...
                'Shepard grid correction step.']);

        if isempty(stnData(:,1))    %If no adjusted or non-adjusted data, print a warning and eliminate correction marker.
            warning('grid_correct:noGHCN', ['No GHCN records '...
                'were found for the downscaling region.  Bias '...
                'correction using Shepard' char(39) 's '...
                'weighting will not occur; however, each iteration ' ...
                'the correction algorithm will look for stations, '...
                'wasting a lot of time.']);
            return
        else
            strGhcnTyp = 'org';  
        end
    end
else
    error('Shep_corr:UnknownArg','Unknown input arguments have been passed to the function.')
end     %End of prelimary GHCN correction load.

if ~isempty(stnData(:,1))   %GHCN stations (either adjusted or non-adjusted) do exist.
    %Create cell-array with same dimensions as GHCN array, but NaN values:
    gridData = init_grid_GHCN(stnData);
else
    return
end

%'GHCN_stn_avail' determines which GHCN stations to use for 
%current time step.
stnUse = GHCN_stn_avail(stnData, sMeta.currTime);
%Format of 'StnUseBool': 
    %1st column is row of each GHCN station in 'stnData'
    %2nd column defaults to zero but if the specific GHCN 
    %station has data for the current year, it takes the value 
    %of that year's row index.   

if ~isempty(stnUse)  %Only continue if station data are available for this time step.   
%     %Write ESRI ASCII format header for subroutine:
%     hdr = ESRI_hdr(sTsScaled.lon, sTsScaled.lat, 'cor');
    %'grid_2_GHCN_form' finds grid values at GHCN station locations 
    %and arranges them in a cell array format that matches 'stnData'.
    gridData = geodata_2_GHCN_form(sTsScaled, gridData, stnData);

    %Use 'gridData' and 'stnData' to create bias correction matrix 
    %using Shepard's weighting:
    [corrMat, corrMeta] = shep_weight(stnData, gridData, sTsScaled.lon, sTsScaled.lat, sMeta.currTime, p, r);

    if corrMeta{1,2} > 0 %Indicates number of stations used more than 0
        dataTsShep = currTs + corrMat;

        %In case bias correction makes pre negative:
        if regexpi(sMeta.currVar,'pre')
            dataTsShep(dataTsShep < 0) = 0;
        end

    else
        warning('grid_correct:noGHCN','No Shepard correction took place because no GHCN stations found for current time step.');
    end
end

dataTsShep = single(dataTsShep);

corrMeta = [corrMeta; {'GHCN data type', strGhcnTyp}];


