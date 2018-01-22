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


%Enter the GHCN ID # of any stations you would like to exclude from the
%analysis:
stnExclude = [];

%If the data files are not in mm or Celsius, this allows for conversion
nScale = [];    %Leave empty if no scaling desired.

%Enter the range of years to be used in calculating stats:
statYrs = [1970,2010];



%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules


%%INNERWORKINGS (DO NOT MODIFY):
uiwait(msgbox(sprintf(['Select the folder with ESRI formatted ASCII ' ...
    'climate data for mean monthly temperature or precipitation. \n']),...
    '(Click OK to Proceed)','modal'));
    dataFolder = uigetdir(pwd,['Select the ESRI '...
        'formatted ASCII climate data.']);
    
metStr = 'The meteorological variable has been autodetected as ';
if ~isempty(regexpi(dataFolder,'pre'))
    metVar = 'pre';
    disp([metStr 'precipitation.']);
elseif ~isempty(regexpi(dataFolder,'tmn')) || ~isempty(regexpi(dataFolder,'tmean'))
    metVar = 'tmn';
    disp([metStr 'mean temperature.']);   
else
    metVar = '';
    iCnt = 0;
    while ~isempty(regexpi(dataFolder,'tmn')) || ~isempty(regexpi(dataFolder,'pre'))
        metVar = input(['Please input the short-hand for the ' ...
            'meteorological variable to be compared: ' char(39) ...
            'pre' char(39) ' or ' char(39) 'tmn' char(39)],'s');
        iCnt = iCnt + 1;
        if iCnt == 5
            disp(['You' char(39) 're stuck in a while loop!']);
        end
    end
end

disp(['You chose ' dataFolder ' as the folder containing ' metVar ...
    ' data.']);

    
%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
sepInd = regexpi(pathScript,filesep);
pathScript = pathScript(1:sepInd(end)-1);
cd(pathScript);
addpath(genpath(pathScript));
    
    
for kk = 1 : 2
    if kk == 1
        cmprType = 'adj';
        strCmpr = 'adjusted';
    elseif kk == 2
        cmprType = 'non';
        strCmpr = 'non-adjusted';
    end

    GHCN_2_grid_cmpr_v2(dataFolder,'N',cmprType,metVar,statYrs,nScale,stnExclude);
        %GHCN_2_grid_cmpr_fun(dirTs{ii},wcStr{ii},cmprType,metVar{ii},statYrs,nScale,stnExclude);

    disp(['Finished processing ' strCmpr ' GHCN stats for ' ...
        dataFolder '.' char(10)]);
end