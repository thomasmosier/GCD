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

clear all

nFolds = 3;



%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules

%%Innerworkings (DO NOT MODIFY):
%Initialize:
foldTs = cell(nFolds,1);
statsOut = cell(nFolds,1);

%Get folder with ts data:
for ii = 1 : nFolds
    if ii == 1
        dirSearch = pwd;
    else
        indSep = regexpi(foldTs{ii-1},filesep);
        dirSearch = foldTs{ii-1}(1:indSep(end-1)-1);
    end
    
    foldTs{ii} = uigetdir(dirSearch, ['Select the directory containing '... 
        'time series data you wish to use for producing a climatology. (' num2str(ii) ' of ' num2str(nFolds) ')']);

    disp([foldTs{ii} ' has been selected for processing.'])
    filesTemp{ii} = dir(fullfile(foldTs{ii},'*.asc'));
    gridFiles{ii} = struct2cell(filesTemp{ii});
    gridFiles{ii} = gridFiles{ii}(1,:);
end

%Initiate waitbar:
hWait = waitbar(0,'Data processing has begun.');

%%DATA PROCESSING
for ii = 1 : nFolds
    %INITIALIZE
    statsOut{ii} = nan(length(gridFiles{ii}),5);
    
    %Create output file:
    fileStat = gridFiles{1}{1};
    indExt = regexpi(fileStat,'_');
    fileStat = [fileStat(1:indExt(end)-1), '_stats.txt'];
    pathStatsOut = fullfile(foldTs{ii},fileStat);
    
    %CREATE HEADER
    if ~isempty(regexpi(foldTs{ii},'tmp')) || ~isempty(regexpi(foldTs{ii},'temp'))
        hdrMet = {'Avg Temp (deg C)', 'Min Temp (deg C)', ...
            'Max Temp (deg C)' ,'Std Temp (deg C)', 'Spatially Integrated Temp (deg C)'};
    elseif ~isempty(regexpi(foldTs{ii},'pre')) || ~isempty(regexpi(foldTs{ii},'pr'))
       hdrMet = {'Avg Pre (mm)', 'Min Pre (mm)', ...
            'Max Pre (mm)' ,'Std Pre (mm)', 'Spatially Integrated Pre (mm)'};
    else
        error('stats:varDetect',[char(39) foldTs{ii} char(39) ' contains an unknown variable.']);
    end
    cellHdr = ['Filename', hdrMet];
    strHdr = blanks(0);
    for kk = 1 : numel(cellHdr)
        if kk ~= numel(cellHdr)
            strHdr = [strHdr char(cellHdr{kk}) ', '];
        else
            strHdr = [strHdr char(cellHdr{kk}) char(10)];
        end
    end

    %Write Header:
    fileID = fopen(pathStatsOut,'w');
    fwrite(fileID, strHdr, 'char');
    
    %Loop over files in current folder
    for jj = 1 : length(gridFiles{ii})
        %LOAD DATA
        [data, hdr, meta] = read_ESRI(fullfile(foldTs{ii},gridFiles{ii}{jj}));
        if jj == 1
            %Calculate area:
            [lat, lon] = ESRI_hdr2geo(hdr,meta);
            area = area_geodata(lon, lat);
        else
            [latTemp, lonTemp] = ESRI_hdr2geo(hdr,meta);
            if ~isequal(lat,latTemp) || ~isequal(lon,lonTemp) 
                error('stats:alignment','Grids do not align.')
            end
        end
        
       %COMPUTE STATS
       %Mean:
       statsOut{ii}(jj,1) = nanmean(nanmean(data.*area)) / nanmean(nanmean(area));
       
       %Min:
       %statsOut{ii}(jj,2) = nanmin(nanmin(data.*area / area));
       [vMn1, iMn1] = nanmin(data.*area);
       [vMn2, iMn2] = nanmin(vMn1);
       statsOut{ii}(jj,2) = vMn2 / area(iMn1(iMn2),iMn2);
       
       %Max:
       [vMx1, iMx1] = nanmax(data.*area);
       [vMx2, iMx2] = nanmax(vMx1);
       statsOut{ii}(jj,3) = vMx2 / area(iMx1(iMx2),iMx2);
       %statsOut{ii}(jj,3) = nanmax(nanmax(data.*area / area));
       
       %Std:
       statsOut{ii}(jj,4) = nanstd(nanstd(data.*area)) / nanmean(nanmean(area));
       
       %Spatially-integrated values:
       statsOut{ii}(jj,5) = nansum(nansum(data.*area)) / nanmean(nanmean(area));
       
       %Convert stats to text and write:
       strStats = [gridFiles{ii}{jj}, ', ' num2str(statsOut{ii}(jj,1)), ', ', ...
           num2str(statsOut{ii}(jj,2)), ', ', num2str(statsOut{ii}(jj,3)), ', ',...
           num2str(statsOut{ii}(jj,4)), ', ', num2str(statsOut{ii}(jj,5)) char(10)];
       fwrite(fileID, strStats, 'char');
       
       %Update waitbar:
        warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
        fracComplt = (ii*(length(gridFiles{ii})-1) + jj)/(length(gridFiles{ii})*nFolds);
        waitbar(fracComplt, hWait, ...
            ['The file ' char(39) char(strrep(gridFiles{ii}{jj},'_','\_')) ...
            char(39) ' is being read and processed.']);
        warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
       
    end
    
    fclose(fileID);
end

delete(hWait);

%%CALCULATE AND WRITE AVERAGE STATS
indSep = regexpi(foldTs{1},filesep);
indUnd = regexpi(gridFiles{1}{1},'_');
indDot = regexpi(gridFiles{1}{1},'\.');
fileStatEns = ['ens' gridFiles{1}{1}(indUnd(1):indDot(end)) 'txt'];
pathStatsEns = fullfile(foldTs{1}(1:indSep(end-1)-1),fileStatEns);
fileIDEns = fopen(pathStatsEns,'w');
fwrite(fileIDEns, strHdr, 'char');
%Find ensemble of stats:
statsEns = zeros(size(statsOut{1}));
for ii = 1 : nFolds
   statsEns = statsEns + statsOut{ii};
end
statsEns = statsEns / nFolds;

%Convert stats to text and write:
for ii = 1 : length(gridFiles{ii})
    strStatsEns = [['ensemble' gridFiles{1}{ii}], ', ' num2str(statsEns(ii,1)), ', ', ...
       num2str(statsEns(ii,2)), ', ', num2str(statsEns(ii,3)), ', ',...
       num2str(statsEns(ii,4)), ', ', num2str(statsEns(ii,5)), char(10)];
    fwrite(fileIDEns, strStatsEns, 'char');
end
fclose(fileIDEns);
