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

function print_geostats_v2(pathGeoStats, stnRec, stnYrs, stnStats)



%'StnYrs' is a vector of the number of years of record at each station.
%'stnRec' is a cell-array containing the stn number ({:,1}), lat ({:,2}),
%lon ({:,3}), elev ({:,4}) for each site.
%'stnStats' is a structure where stnStats.data is a matrix of station
%statistics and stn.name is an array of strings which describe the
%statistics used.

nPrec = 4;

fidGeoStats = fopen(pathGeoStats, 'w');

%%PRINT HEADER:
strGeoHdr = {'ghcn_id,lat,lon,elev,n_years'};

for ii = 1 : length(stnStats)
    statLab = stnStats(ii).name;
    spc2Undr = strfind(statLab, ' ');
    if ~isempty(spc2Undr)
        statLab(spc2Undr) = '_';
    end
    
    strGeoHdr = [char(strGeoHdr), ',', statLab];
end

fprintf(fidGeoStats, '%s\n', char(strGeoHdr));

%%WRITE STATION STATS:
%This creates a string which is used in fprintf to indicate data formatting.
strDatFlt = ['%.', num2str(nPrec), 'f'];
strDatInt = '%u';
strDataWrt = strDatInt;
for ii = 1 : length(stnStats) + 4   %Iterate over hdr items and statistics
    
    if ii == 4  %n_years should be integer
        strDataWrt = [char(strDataWrt), ',', char(strDatInt)];
    else
        strDataWrt = [char(strDataWrt), ',', char(strDatFlt)];
    end
end

%This isn't necessary anymore:
%strDataWrt = [char(strDataWrt), '\n' char(10)];

%For testing:
%disp(['The GeoStats string being written is:' strDataWrt])

for ii = 1 : length(stnYrs(:))
    geoStatsVec = [stnRec{ii,1}, stnRec{ii,2}, stnRec{ii,3}, stnRec{ii,4}, ...
        stnYrs(ii)];
    for jj = 1 : length(stnStats)
        geoStatsVec = [geoStatsVec, stnStats(jj).data(ii)];
    end
    
    fprintf(fidGeoStats, strDataWrt, geoStatsVec);
    fprintf(fidGeoStats,'\n');
end

fclose(fidGeoStats);

end