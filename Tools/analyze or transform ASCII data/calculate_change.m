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



%This script calculates the percent and absolute change between two ESRI ASCII files


uiwait(msgbox(sprintf(['You will be prompted to select two files.  '...
    'This script will calculate the percent change between these '...
    'two files.  Both files must have the same spatial extent.']),...
    '(Click OK to Proceed)','modal'));

[fileFin, pathFin] = uigetfile({'*.asc';'*.txt'}, ...
    'Select the ESRI ASCII file for the final state.');
fileFin = fullfile(pathFin,fileFin);
disp([char(39) fileFin char(39) ' has been selected as the final state.']);

indSep = regexpi(pathFin, filesep);
pathSrch = pathFin(1:indSep(end-1)-1);
[fileDebut, pathDebut] = uigetfile({'*.asc';'*.txt'}, ...
    'Select the ESRI ASCII file for the final state.', pathSrch);
fileDebut = fullfile(pathDebut,fileDebut);
disp([char(39) fileDebut char(39) ' has been selected as the initial state.']);

[dataFin, hdrFin, ~] = read_ESRI(fileFin);
[dataDebut, hdrDebut, ~] = read_ESRI(fileDebut);

if ~isequal(hdrFin,hdrDebut)
    error('percent_change:misalign','The two files do not align.');
end

dataPerChng = 100*(dataFin - dataDebut)./dataDebut;
dataPerChng(dataDebut == 0) = NaN;

dataAbsChng = dataFin - dataDebut;
dataAbsChng(dataDebut == 0) = NaN;

indExt = regexpi(fileFin,'\.');
filePerChng = [fileFin(1:indExt(end)-1), '_relative_change.asc'];
disp(ln80(['The percent change data are being written to ' char(39) ...
    filePerChng char(39) '.']));
fileAbsChng = [fileFin(1:indExt(end)-1), '_absolute_change.asc'];
disp(ln80(['The absolute change data are being written to ' char(39) ...
    fileAbsChng char(39) '.']));

%Write data:
write_ESRI_v3(dataPerChng, hdrFin, filePerChng, 1);
write_ESRI_v3(dataAbsChng, hdrFin, fileAbsChng, 1);
    