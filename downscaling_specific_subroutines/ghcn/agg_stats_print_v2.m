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


function agg_stats_print_v2(fidStats, aggStats, statSpec, nStn)



strAggHdr = 'The statistics, aggregated over all parameters, are:';
fprintf(fidStats, '%s', strAggHdr);

%'length(stnStats)' is used here to record the total number of contributing
%stations, but because the stats are aggregated, the length of the entry is
%one (i.e only one lien of stats will be written).
aggHead = [num2str(nStn), ' Stns Used'];


strHdr = GHCN_stat_hdr(aggHead, 'Aggregated', aggStats, statSpec);
% disp(repMnthHdrStr);
fprintf(fidStats, '%s', [char(10) strHdr]);

% nStats = length(stnStats);
% 
% hdrData = cell(nStats, 1);
% for ii = 1 : nStats
%    hdrData{ii} = char(stnStats(ii).name);
% end


GHCN_stat_data(aggHead, aggStats, statSpec, fidStats);

end