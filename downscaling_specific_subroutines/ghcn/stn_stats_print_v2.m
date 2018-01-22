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

function stn_stats_print_v2(fidReport, ghcnName, stnStats, statSpec)



%DESCRIPTION:
repStnInfoStr = ['The following statistics are organized by GHCN station ID ' char(10) ...
    '(the accompanying CSV .txt file provides these stats in a format importable in ArcMap).'];
fprintf(fidReport, '%s', [repStnInfoStr char(10)]);


%%PRINT HEADER:
if isempty(ghcnName)
    strHdr = 'No GHCN station data available for this combination of time elements and region';

    fprintf(fidReport, '%s', char(strHdr));
else
    strHdr = GHCN_stat_hdr(ghcnName{1},'GHCN ID',stnStats, statSpec);

    fprintf(fidReport, '%s', char(strHdr));

    %%PRINT STATS:
    hdr = [ghcnName(:); 'Avg'];
    GHCN_stat_data(hdr, stnStats, statSpec, fidReport);
end
end