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

function write_GHCN(data, pathData)



[fileRoot, ~, ~ ] = fileparts(pathData);
if ~exist(fileRoot,'dir')
   mkdir(fileRoot); 
end

if isunix
    pathData = [filesep pathData];
end

hdrLbl = 'GHCN ID, Lat, Lon, elevation';
hdrMnths = 'Year, Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec';
fiddata = fopen(pathData, 'w');

for ii = 1 : length(data(:,1))
    fprintf(fiddata,'%s \r\n', hdrLbl);
    fprintf(fiddata,'%i, %.4f, %.4f, %i\r\n',  data{ii,1:4});
    
    fprintf(fiddata,'%s\r\n', hdrMnths);
    for jj = 1 : length(data{ii,5}(:,1))
        fprintf(fiddata,'%i, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f \r\n', data{ii,5}(jj,:));
    end
    
    fprintf(fiddata, '\r\n');
    
end

fclose(fiddata);

end