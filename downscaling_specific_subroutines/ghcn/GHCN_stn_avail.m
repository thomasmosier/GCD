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

function stnUseBool = GHCN_stn_avail(stnData, gridDate)


%Format of 'StnUseBool': 
    %1st column is row of each GHCN station in 'stnData'
    %2nd column defaults to zero but if the specific GHCN station has
    %data for the current year, it takes the value of that year's row
    %index.    

stnUseBool = (1:length(stnData(:,1)));

for jj = length(stnData(:,1)) : -1 : 1
    %Determine if each station has data for current year-month and, if yes,
    %save row index.
    stnRowTemp = find(stnData{jj,5}(:,1) == gridDate(1));
    if ~isempty(stnRowTemp) && ~isnan( stnData{jj,5}(stnRowTemp, 1 + gridDate(2)) )
        continue
    else
        stnUseBool(jj) = [];
    end
end