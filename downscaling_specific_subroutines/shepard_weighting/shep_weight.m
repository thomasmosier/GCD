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

function [corrMat, corrMeta] = shep_weight(stnData, gridData, lon, lat, date, p, r)



%This function uses Shepard's weighting (see 1968 paper) to produce a
%correction matrix at each time-step, where the list of correction factors
%is the difference between GHCN station records and downscaled data
%surface.

%CorrMeta is equal to number of stations used in Shepard's weighting 
%(=0 is no correction occured due to lack of stations.)

%'stnData' and 'gridData' are cell arrays with the form:
    %{i,1} is the GHCN station ID,
    %{i,2} is the cell-center latitude, 
    %{i,3} is the cell-center longitude,
    %{i,4} is the station elevation, and 
    %{i,5} is an m by 13 cell-array with yearly data (first column is year).

%'ep' provides significant figures for rounding.
ep = 1000;

%Initiliaze CorrMeta:
corrMeta = cell(3,2);
    corrMeta{1,1} = 'Number of Stations Used';
    corrMeta{2,1} = 'Mean of Gridded Bias';
    corrMeta{3,1} = 'Standard Deviation of Gridded Bias';

if date(1) > date(2)
    mnth = date(2);
    yr = date(1);
elseif date(1) < date(2)
    mnth = date(1);
    yr = date(2);
else
    error('shep_weight:incorrectDates',['The elements of the date vector are ' num2str(date(1)) ...
        ' and ' num2str(date(2)), ' which cannot be correct.']);
end
    
[lonGrid, latGrid] = meshgrid(lon, lat);

lonGrid = round(lonGrid*ep)/ep;
latGrid = round(latGrid*ep)/ep;


%Initialize:
gridUse = nan(length(gridData(:,1)),3);
stnUse  = nan(length(gridData(:,1)),3);


%Find stations with meteorological values for the current time step:
useCntr = 1;    %Cntr for loop over meteorological stations
for ii = 1 : length(gridData(:,1))
    currYr = find( gridData{ii,5}(:,1) == yr );
    if ~isempty(currYr) && ~isnan(gridData{ii,5}(currYr, mnth+1))
        %Copy latitudes:
        gridUse(useCntr,1) = gridData{ii,2};
         stnUse(useCntr,1) =  stnData{ii,2};

        %Copy longitudes:
        gridUse(useCntr,2) = gridData{ii,3};
         stnUse(useCntr,2) =  stnData{ii,3};  

        %Copy meteorological values:
        gridUse(useCntr,3) = gridData{ii,5}(currYr, mnth+1);
         stnUse(useCntr,3) =  stnData{ii,5}(currYr, mnth+1);

        useCntr = useCntr + 1;
    end
end

%delete empty rows for initialized matrices:
if useCntr <= length(stnUse(:,1))
    gridUse(useCntr:end,:) = [];
     stnUse(useCntr:end,:) = [];
end

%Initialize correction matrix
corrMat = nan( length(lat), length(lon) );

if ~isempty( stnUse(:,1) )    %Only proceed if at least one station exists for current timestep.
    for kk = 1 : length( corrMat(:,1) ) %interate over rows of corrMat
        for ll = 1 : length( corrMat(1,:) ) %iterate over columns of corrMat

            %Check if there is a station within the current grid cell 
            stnInCell = intersect( ...
                find( round(gridUse(:,1)*ep)/ep == latGrid(kk,1) ), ...
                find( round(gridUse(:,2)*ep)/ep == lonGrid(1,ll) ) ...
                );

            if ~isempty(stnInCell)  %If grid pt contains stn, correction factor is station bias: 
                corrMat(kk,ll) = stnUse(stnInCell, 3) - gridUse(stnInCell, 3);

            else %If grid pt does not contain a station, compute Shepard's weight:
                deltaLat = stnUse(:,1) - latGrid(kk,1);   %Vector
                deltaLon = stnUse(:,2) - lonGrid(1,ll);   %Vector

                %Vector of distances between current cell and each
                %station with value at current time step.
                dJK = ...
                    r*2*asin(sqrt( sind(deltaLat/2).^2 ...
                    + cosd(latGrid(kk,1)).*cosd(stnUse(:,1)).*sind(deltaLon/2).^2 ));

                corrMat(kk,ll) ...
                    = sum( (dJK.^(-p)).*(stnUse(:,3) - gridUse(:,3)) )...
                    /sum(dJK.^(-p));
            end

        end
    end

    corrMeta{1,2} = length(stnUse(:,1)); %Number stations used
    corrMeta{2,2} = mean(mean(corrMat)); %Mean bias over region
    corrMeta{3,2} = std(std(corrMat)); %Standard deviation of bias  
else    %Case where no stn data exists for current time step.
    corrMeta{1,2} = 0; %Number stations used
    corrMeta{2,2} = nan; %Mean bias over region
    corrMeta{3,2} = nan; %Standard deviation of bias
end

end     %End of function!