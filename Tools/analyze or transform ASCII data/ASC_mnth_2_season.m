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

%%User settings:
nFold = 2; %Can input multiple folders to be processed in series  

mnths = {(1:12), [3,4,5], [6,7,8,9], [10,11,12,1,2]};  %Creates monthly climatologies for each month in 'mnths' vector
calcTyp = 'mean'; %Either 'sum' (sums data within month vector) or 'mean' (averages data within month vector).
    sclTs = 1;  %Should be 1 except for PRISM temperature
%     calcAnom = 0;
    
    
%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules


%%Innerworkings (DO NOT MODIFY):
% if calcAnom == 1
%     warning('season:noAnom','This function currently does not have teh ability to calculate seasonal anomalies.')
%     disp(['The anomaly calculation option has been selected.  '...
%         'The user will be required to choose two folders.  ' ...
%         'The first contains the monthly time-series or climatologies '...
%         'and the second contains the reference data against which ' ...
%         'the anomaly will be calcualted.']);
% end

%Get folder with ts data:
foldTs = cell(nFold,1);
for kk = 1 : nFold
    if kk == 1
        pathStart = pwd;
    else
        indSep = regexpi(foldTs{kk-1},filesep);
        pathStart = foldTs{kk-1}(1:indSep(end-1)-1);
    end

    foldTs{kk} = uigetdir(pathStart,['Select the directory containing '... 
        'ESRI gridded ASCII monthly data (' num2str(kk) ' of ' num2str(nFold) ')']);
    disp([foldTs{kk} ' has been selected as the folder of monthly ESRI ASCII files.']);
end
disp(blanks(0));

% %Get folder with anom data:
% if calcAnom == 1
%     uiwait(msgbox(sprintf(['Select the directory containing '... 
%         'ESRI gridded ASCII baseline data to use for calculating ' ...
%         'the anomaly. \n']), '(Click OK to Proceed)','modal'));
%     foldAnom = uigetdir(['Select the directory containing '... 
%         'ESRI gridded ASCII baseline data to use for calculating the '...
%         'anomaly.']);
% 
%     disp([foldAnom ' has been selected as the folder of monthly ESRI ASCII files for the reference.']);
% 
%     filesAnomTemp = dir(fullfile(foldAnom,'*.asc'));
%     gridAnomFiles = struct2cell(filesAnomTemp);
%     gridAnomFiles = gridAnomFiles(1,:);
% end

for kk = 1 : numel(foldTs(:))
    for ii = 1 : numel(mnths(:))
        sMeta.mnths = mnths{ii};
        [foldMainOut, foldStdOut] = season(foldTs{kk},sMeta,calcTyp,sclTs);
    end
    disp('');
end