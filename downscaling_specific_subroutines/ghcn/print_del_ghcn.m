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

function print_del_ghcn(fidReport, recStnNan, reason)



if ~isempty(recStnNan)
    if regexpi(reason,'nan')
        stnNanHdrStr = ['The following stations were excluded from analysis ' ... 
            'because the gridded data was entirely NaN at the location.'];
    elseif regexpi(reason,'years')
        stnNanHdrStr = ['The following stations were excluded from analysis ' ... 
        'because they did not contain a sufficient number of years of data.'];
    else
        stnNanHdrStr = ['The following stations were excluded from analysis ' ... 
        'for an unundefined reason.'];
    end
%     disp(stnNanHdrStr);
    fprintf(fidReport, '%s', [stnNanHdrStr char(10)]);
    
    for ii = 1 : length(recStnNan)
        stnNaNstr = [num2str(recStnNan{ii}) ', '];
        fprintf(fidReport, '%s', [stnNaNstr char(10)]);
%         disp(stnNaNstr);
    end
    fprintf(fidReport, '%s', char(10));
end

end