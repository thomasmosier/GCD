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

nFold = 2;


%Iniatiate variables:
fileESRI = cell(nFold,1);
foldESRI = fileESRI;
pathESRI = fileESRI;

%GUI to locate files:
for jj = 1 : nFold
    if jj == 1
        pathStrt = pwd;
    else
        indSep = regexpi(foldESRI{jj-1},filesep);
        pathStrt = foldESRI{jj-1}(1:indSep(end));
    end
    [fileESRI{jj}, foldESRI{jj}, ~] = uigetfile({'*.asc'}, ...
        ['Select ESRI ASCII file ' num2str(jj) ' of ' num2str(nFold) '.'], pathStrt);
    pathESRI{jj} = fullfile(foldESRI{jj}, fileESRI{jj});
end

%Read data and calculate average:
for jj = 1 : nFold
   [data,hdr{jj},meta{jj}] = read_ESRI(pathESRI{jj});
    
   if jj == 1
       dataAvg = zeros(size(data));
   elseif jj ~= 1 && ~isequal(hdr{jj},hdr{jj-1})
        error('ESRI_avg:hdrmisalign',['Hdr' num2str(jj) ' and ' num2str(jj-1) ' are not equal.']);
   end
    
   dataAvg = dataAvg + data;
end

dataAvg = dataAvg / nFold;

%Write file:
indExt = regexpi(pathESRI{1}, '\.');

pathOut = [pathESRI{1}(1:indExt(end)-1), '_ESRIavg.asc' ];
write_ESRI_v4(dataAvg,hdr{1},pathOut,0);
disp(['Average of ESRI ASCII files written to ' pathOut]);
