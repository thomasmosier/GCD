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


%Read NetCDF files, produce climatologies, and write to ASCII

%%ENTER INFO:
%Choose to create either time-series or climatologies:
wrtType = 'ts'; %Either 'ts' or 'clim' depending on which is needed.

%Enter metadata info:
% sMeta.currVar = 'pre'; %'pre' 'tmn' 'tmp', 'tmx'
yrs = [2020, 2020]; %Select range of years
sMeta.mnths = (1:1); %Selects months
sMeta.crd = [-180, 180, -89, 89]; %Select bounding box


%%DO NOT EDIT:
%Other metadata:
sMeta.hisTS = '';
if ~isempty(regexpi(wrtType,'c'))
    sMeta.yrsOut = [NaN, NaN];
    sMeta.yrsClim = yrs; 
elseif ~isempty(regexpi(wrtType,'ts'))
    sMeta.yrsOut = yrs;
    sMeta.yrsClim = [NaN, NaN]; 
end

%Get the NetCDF file to be converted to ASCII:
% uiwait(msgbox(sprintf(['Select the NetCDF file for ' sMeta.currVar ' to be converted to '...
%     'ASCII \n']), '(Click OK to Proceed)','modal'));
[fileData, dirData, ~] = uigetfile({'*.nc';'*.grb'}, 'Select the NetCDF file');

pathData = fullfile(dirData,fileData);

if ~isempty(regexpi(fileData,'tas')) || ~isempty(regexpi(fileData,'ts')) 
    sMeta.currVar = 'tmp';
elseif ~isempty(regexpi(fileData,'pre')) || ~isempty(regexpi(fileData,'pr'))
    sMeta.currVar = 'pre';
else
    error('NC_2_ASC',ln80(['The parameter for ' dataFile ...
        ' is unknown.  Add code for this case.']));
end
%If this is the first time the script has been run, expand search path and
%load NC Toolbox.
if ~exist('ncLd','var') 
    %Add path of modules that may be called:
    pathScript = mfilename('fullpath');
    [pathScript, ~, ~] = fileparts(pathScript);
    indParent = regexpi(pathScript, filesep);
    indParent = indParent(end-1); %Assumes scripts are in distribution hierarchy.
    addpath(genpath(pathScript(1:indParent-1))); %Add modules
    
    %Find NC Toolbox:
    ncSearch = genpath(pathScript(1:indParent-1));
%     ncSearch = path;
    indNC = regexpi(ncSearch,'nctoolbox-nctoolbox');
    indLine = regexpi(ncSearch,';');
    indNCUse(1) = find(indLine < indNC(1),1,'last');
    indNCUse(2) = find(indLine > indNC(1),1,'first');
    %Run NC Toolbox setup:
    run(fullfile(ncSearch(indLine(indNCUse(1))+1 : indLine(indNCUse(2))-1),'setup_nctoolbox.m'));
    
    %Record that setup has occured:
    ncLd = 1;
end

%Initiate waitbar:
hWait = waitbar(0,'Data processing has begun.');

for ii = 1 : length(sMeta.mnths)
    sMeta.currTime = [NaN,sMeta.mnths(ii)];

    sData = read_geodata(pathData, sMeta, wrtType);

    %Ensure sData.data has 3 indices:
    if length(size(sData.data)) == 2
        sData.data = reshape(sData.data, [1, size(sData.data)]);
    end
    
    %Check that GCM grid cells conform to specified bounds.  If not,
    %resample:
    if ii == 1
        up = 0;
        
        if round(10*(sMeta.crd(1) - sData.lon(1)))/10 > 0
            up = 1;
        elseif round(10*(sMeta.crd(2) - sData.lon(end)))/10 > 0
            up = 1;
        elseif round(10*(sMeta.crd(3) - sData.lat(end)))/10 > 0
            up = 1;
        elseif round(10*(sMeta.crd(4) - sData.lat(1)))/10 > 0
            up = 1;
        end
        
        up = 0;
%         if up == 1
%             strUp = input(['The input data and the bounding ' ...
%                 'coordinates specified do not align.  ' char(10) ...
%                 'Would you like to resample the input data so that ' ...
%                 'they do align? (' char(39) 'Y' char(39) ' or ' ...
%                 char(39) 'N' char(39) ')' char(10)],'s');
%             if ~isempty(regexpi(strUp,'y'))
%                 up = 1;
%             else
%                 up = 0;
%             end
%         end
    end
    
    if up == 1
        if ii == 1
            %Find even grid spacing for lat-lon that conforms to specificed
            %input region
            dLon = abs(sMeta.crd(2)-sMeta.crd(1));
            dLat = abs(sMeta.crd(4)-sMeta.crd(3));
            
            divisor = 1:100000;
            divLon = dLon ./ divisor;
            divLat = dLat ./ divisor;
            
            threshold = 0.005; %Threshold for saying that grid spacing is close enough;
            
            remStep = rem(divLat,divLon); 
            nStep = find(remStep < threshold,1,'first');
            
            sRef.lon = linspace(sMeta.crd(1),sMeta.crd(2),nStep);
            sRef.lat = linspace(sMeta.crd(4),sMeta.crd(3),nStep)';
            
            stpLatOrg = mean(abs(diff(sData.lat)));
            stpLonOrg = mean(abs(diff(sData.lon)));
            stpLatNew = mean(abs(diff(sRef.lat)));
            warning('NC_2_ASCI:resample',['The NetCDF data is beign '...
                'resampled from an original spatial resolution of ' ...
                char(10) num2str(stpLatOrg) ' (lat) and ' ...
                num2str(stpLonOrg) ' (lon) to' char(10) ...
                num2str(stpLatNew) ' (uniform in lat and lon)' ...
                char(10) ' because of the requirements that the grid '...
                'resolution be uniform in lat and lon ' char(10) ...
                'and that the exterior bounds conforms to the specified inputs.']);
        end
        sData = regrid_geodata(sData,sRef);
        sData.data = sData.dataRe;
        sData.lat = sData.latRe;
        sData.lon = sData.lonRe;
    end
    
    %Write each grid to ASCII file:
    %Find NetCDF name root to use in output naming:
    indU = regexpi(fileData,'_');
    nmData = fileData(1:indU(end));
    
    %Create output name:
    if ~isempty(regexpi(wrtType,'c'))
        pathWrite = fullfile(dirData, nmData(1:end-1), ...
            [nmData, num2str(sMeta.yrsClim(1)), 'thru', ...
            num2str(sMeta.yrsClim(2)), '_' num2str(sMeta.mnths(ii)), '.asc']);
            %Write data:
        write_geodata(pathWrite,sData,sMeta,1,'ascii');
    elseif ~isempty(regexpi(wrtType,'ts'))
        for jj = 1 : (yrs(2) - yrs(1) + 1)
            sMeta.currTime = [(sMeta.yrsOut(1)+jj-1), sMeta.mnths(ii)];
            pathWrite = fullfile(dirData, nmData(1:end-1), ...
                [nmData, num2str(sMeta.currTime(1)), ...
                '_' num2str(sMeta.mnths(ii)), '.asc']);
            %Write data:
            write_geodata(pathWrite,sData,sMeta,1,'ascii');
        end
    end



    %Update waitbar:
    warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    fracComplt = ii/length(sMeta.mnths);
    waitbar(fracComplt, hWait, ...
        ['The file ' char(39) char(strrep(pathWrite,'_','\_')) ...
        char(39) ' is being written.']);
    warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
end

delete(hWait);