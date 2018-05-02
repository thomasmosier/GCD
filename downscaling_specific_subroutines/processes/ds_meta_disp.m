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

function ds_meta_disp(sPath, dataLd, indSim, sMeta, varCurr)


%Display version of current variable, with first letter capitalized
metVarCap = [upper(varCurr(1)), varCurr(2:end)];

%DISPLAY USER CHOICES
disp([metVarCap ' data for the years ' num2str(sMeta.yrsOut(1)) ...
    ' thru ' num2str(sMeta.yrsOut(2)) ' and months ' ...
    num2str(sMeta.mnthsOut(1)) ' to ' num2str(sMeta.mnthsOut(end))  ...
    ' will be downscaled.']);

disp(['The baseline years (e.g. those used in the historic reference climatology) are ' ...
    num2str(sMeta.yrsBase(1)) ' thru ' num2str(sMeta.yrsBase(2)) '.']);

fldsDisp = fieldnames(sPath);

for ii = numel(fldsDisp) : -1 : 1
    if strcmpi(fldsDisp{ii}(1:2), 'nm')
        fldsDisp(ii) = [];
    end
end
clear ii

fldsDisp(strcmpi(fldsDisp, 'output')) = [];

for ii = 1 : numel(fldsDisp)
    pathCurr = sPath.(fldsDisp{ii});
    if iscell(pathCurr)
        if numel(pathCurr(:)) == 1
            pathCurr = char(pathCurr);
        else
            pathCurr = char(pathCurr{indSim});
        end
    end
    disp([dataLd{ii,2} ' ' dataLd{ii,1} ' data are being loaded from ' pathCurr '.']);
end
clear ii
    
disp(['The ' char(sMeta.methodDisp) ' method will be implemented, '...
    'where the low-resolution grids are being interpolated with ' ...
    char(sMeta.intrp) '.' char(10)]);
if regexpbl(sMeta.method,'his') && sMeta.shepard == 1
    disp(['The downscaled data will be bias corrected against available station'...
        ' data from the Global Historical Climatology Network, '...
        'resulting in an additional downscaled time-series output.']);
end