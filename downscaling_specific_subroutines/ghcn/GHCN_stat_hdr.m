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

function strHdr = GHCN_stat_hdr(strCateg, hdrStrt, stnStats, statSpec)



%This function writes the header string for any of the GHCN statistics
%calculated.  Main use is inside '*_stats_print_v2.m'

%'hdrStrt' is the name of the category (i.e. 'Month', 'Agreggated', etc.)
%'stnStats' is a structure with the stat title under '.name' and the
%statistics under '.data'.
%'statSpec' indicates the data writing format (e.g. '%09.5f').

extrSpc = '';

    warning('OFF','MATLAB:nonIntegerTruncatedInConversionToChar');
nAddSpc = 2 + length(char(strCateg)) - length(char(hdrStrt));
    warning('ON','MATLAB:nonIntegerTruncatedInConversionToChar');
if nAddSpc > 0
    for ii = 1 : nAddSpc
        extrSpc = [extrSpc, ' '];
    end
end

strHdr = [extrSpc, hdrStrt];

nStats = length(stnStats);

%Find number of spaces needed between header entries (based on output
%formatting):
indFrm = strfind(statSpec, char(46));
nHdrSpc = str2double(statSpec(2:indFrm-1));

for ii = 1 : nStats
    strHdr = [char(strHdr) ', '];
    nStatSpc = length(char(stnStats(ii).name));
    for jj = 1 : nHdrSpc - nStatSpc
        strHdr = [char(strHdr), char(32)];
    end
    strHdr = [char(strHdr), char(stnStats(ii).name)];
end