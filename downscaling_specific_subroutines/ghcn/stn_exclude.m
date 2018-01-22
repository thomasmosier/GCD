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

function [recStnEx, gridData, stnData, stnYears] = stn_exclude(gridData, stnData, stnYears, stnExclude)



%Uses GHCN station ID numbers in 'stnExclude' to find and remove unwanted
%station data points from 'gridData', 'stnData', and 'stnYears'.

recStnEx = [];
if ~isempty(stnExclude)
    for ii = 1 : sum(~isnan(stnExclude))
        stnVec = cell2mat(gridData(:,1));
        stnDel = find(stnVec == stnExclude(ii));
        if ~isempty(stnDel)
            recStnEx = cat(1, recStnEx, gridData(stnDel,1));
            gridData(stnDel,:) = [];
            stnData( stnDel,:) = [];
            stnYears(stnDel,:) = [];
        end
    end
end

end