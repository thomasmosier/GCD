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

function gridArray = geodata_2_GHCN_form(sGeodata, gridArray, ghcnArray, varargin)

varUse = 'data';
if ~isempty(varargin(:))
   for ii = 1 : numel(varargin(:))
       if strcmpi(varargin{ii}, 'var')
           varUse = varargin{ii+1};
       end
   end
end

[dateRef, gcmUnits] = NC_time_units(sGeodata.attTime);
if ~regexpbl(gcmUnits, {'days','since'})
    error('geodata_2_GHCN_form:timeUnits',['The time units of the '...
    	'current gridded data are ' gcmUnits ', which are not recognized.']) 
end
     
datesGrid = days_2_date(sGeodata.time, dateRef, NC_cal(sGeodata.attTime));
 
if isfield(sGeodata,'lat')
    varLat = 'lat';
elseif isfield(sGeodata,'latitude')
    varLat = 'latitude';
end
latEdg = box_edg(sGeodata.(varLat));
latEdg = [latEdg(1), latEdg(end)];

if isfield(sGeodata,'lon')
    varLon = 'lon';
elseif isfield(sGeodata, 'longitude')
    varLon = 'longitude';
end

lonEdg = box_edg(sGeodata.(varLon));
lonEdg = [lonEdg(1), lonEdg(end)];


for ii = 1 : numel(sGeodata.time)
    stnUse = GHCN_stn_avail(ghcnArray, datesGrid(ii,:));
    
    if isempty(stnUse)
       continue 
    end

    %This for loop iterates over all stations present for a single
    %meteorological grid
    for jj = 1 : numel( stnUse )
        if ghcnArray{stnUse(jj),3} < lonEdg(1) || ghcnArray{stnUse(jj),3} > lonEdg(end)
            warning('geodata_2_GHCN_form:lonOut',['The longitude of '...
                'the current GHCN station falls outside the gridded data.']);
            continue
        end
        if ghcnArray{stnUse(jj),2} > latEdg(1) || ghcnArray{stnUse(jj),2} < latEdg(end)
            warning('geodata_2_GHCN_form:latOut',['The latitude of the '...
                'current GHCN station falls outside the gridded data.']);
            continue
        end
        
            
        [~, indLat] = min(abs(ghcnArray{stnUse(jj),2} - sGeodata.(varLat)));
        [~, indLon] = min(abs(ghcnArray{stnUse(jj),3} - sGeodata.(varLon)));

        %gridData is a cell array where: 
        %{i,1} is the GHCN station ID,
        %{i,2} is the cell-center latitude, 
        %{i,3} is the cell-center longitude,
        %{i,4} is the station elevation, and 
        %{i,5} is an m by 13 cell-array with yearly data (first column is year).

        %Populate grid lat/lon
            %?Add elevation functionality later
        gridArray{stnUse(jj), 2} = sGeodata.(varLat)(indLat);
        gridArray{stnUse(jj), 3} = sGeodata.(varLon)(indLon);

        if numel(size(sGeodata.(varUse))) == 2
            temp = sGeodata.(varUse);
            sGeodata.(varUse) = nan([1,size(temp)]);
            sGeodata.(varUse)(1,:,:) = temp;
        end
        indYr = find(datesGrid(ii,1) == gridArray{stnUse(jj), 5}(:,1));
        if ~isempty(indYr)
            gridArray{stnUse(jj), 5}(indYr, 1+datesGrid(ii,2)) = sGeodata.(varUse)(ii, indLat, indLon);
        end
    end


end