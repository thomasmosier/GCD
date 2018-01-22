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


%Designed to be used for converting time series downscaled output from
%associated downscaling script to climatologies, as specified in 'user
%settings' section.

clear all
clc

nFold = 2;
%%User settings:
climYrs = {[1976,2005], [1979,2008], [1982,2011]};  %Defines years to be used in creating climatology.
mnths = (1:12);  %Creates monthly climatologies for each month in 'mnths' vector
calcTyp = 'mean'; %'mean' or 'sum'
scl = 1;  %Should be 1 except for PRISM temperature


%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules

%%Innerworkings (DO NOT MODIFY):
%Get folder with ts data:
foldTs = cell(nFold,1);
for jj = 1 : nFold
    if jj == 1
       startpath = pwd; 
    else
        indStart = regexpi(foldTs{jj-1},filesep);
        startpath = foldTs{jj-1}(1:indStart(end-1)-1);
    end
    foldTs{jj} = uigetdir(startpath,['Select the directory of time-series data to ' ...
        'produce climatology from (' num2str(jj) ' of ' num2str(nFold) ')']);
    disp([foldTs{jj} ' has been selected as time-series folder ' num2str(jj) ' of ' num2str(nFold) '.'])
end
disp(char(10));


if iscell(climYrs)
    nClim =  numel(climYrs(:));
else
    nClim = 1;
end

% hWait = waitbar(0,'Climatology processing has begun.');
for kk = 1 : nClim
    for jj = 1 : nFold
        if iscell(climYrs)
            currYrs =  climYrs{kk};
        else
            currYrs = climYrs;
        end
        
        [foldMainOut, foldStdOut] = climatology(foldTs{jj}, currYrs, mnths, calcTyp, scl);
    end
end
% delete(hWait);