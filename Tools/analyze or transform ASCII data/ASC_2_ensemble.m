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


%Create ensemble average and standard deviation of geodata in same
%reference frame.
nEns = 8;
nFiles = 6;

%This parameter defines the maximum resolution if grids are resampled to
%common dataframe:
maxRes = 0.0083; %Approximately 30-arcseconds

%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules


fileUse = cell(nEns,nFiles);
path = cell(nEns,nFiles);
dirOut = cell(nEns,1);

uiwait(msgbox(sprintf(['You will be prompted to select ' ...
    num2str(nFiles) ' files.  If the selected files do not have the '...
    'same extent an option will be presented to find a common grid.']),...
    '(Click OK to Proceed)','modal'));
%Loop to get all file paths:
for kk = 1 : nEns
    for ii = 1 : nFiles
        if kk == 1
            if ii == 1
               pathSrch = pwd;
            else
                indSep = find(path{kk,ii-1},filesep);
                pathSrch = path{kk,ii-1}(1:indSep(end-1)-1);
            end
        else 
            indSep = find(path{kk-1,ii},filesep);
            pathSrch = path{kk-1,ii}(1:indSep(end)-1);
        end
        [fileUse{kk,ii}, path{kk,ii}] = uigetfile({'*.asc';'*.txt'},['Select the ESRI ASCII file for ensemble ' num2str(kk) ' (file ' ...
            num2str(ii) ' of ' num2str(nFiles) ')'],pathSrch);
        fileUse{kk,ii} = fullfile(path{kk,ii},fileUse{kk,ii});
        disp([char(39) fileUse{kk,ii} char(39) ' has been selected as file ' num2str(ii) '.']);
    end
    disp('');
end

%%Loop to process all files:
%Initialize arrays
data = cell(size(fileUse));
hdr = cell(size(fileUse));
meta = cell(size(fileUse));
lat = cell(size(fileUse));
lon = cell(size(fileUse));

%Load data
for kk = 1 : nEns
    up = 0;
    for ii = 1 : nFiles
        [data{kk,ii}, hdr{kk,ii}, meta{kk,ii}] = read_ESRI(fileUse{kk,ii});

        %Check current file aligns with previous:
        if ii~= 1
            if ~isequal(hdr{kk,ii},hdr{kk,ii-1}) && up == 0;
                strIn = input([char(39) fileUse{kk,ii-1} char(39) ' and ' ...
                    char(39) fileUse{kk,ii} char(39) ' do not align.' ...
                    'Would you like to interpolate these files to a common grid? (Y or N)'],'s');
                if ~isempty(regexpi(strIn,'y'))
                    up = 1;
                else
                    error('ASC2Ensemble:misalign',['Two of the files do '...
                        'not align and a common dataframe has not been chosen.']);
                end
            end
        end 
    end

    %If data have different spatial properties, create common grid
    if up == 1
        latAll = [];
        lonAll = [];
        %Load all input grids:
        for ii = 1 : nFiles
            [lat{kk,ii}, lon{kk,ii}] = ESRI_hdr2geo(hdr{kk,ii},meta{kk,ii});
            latAll = [latAll; reshape(lat{kk,ii},[],1)];
            lonAll = [lonAll, reshape(lon{kk,ii},1,[])];
        end
        latAll = sort(latAll,'descend');
        lonAll = sort(lonAll,'ascend');
        dLatAll = diff(latAll);
        dLonAll = diff(lonAll);
        minDLat = min(abs(dLatAll));
        minDLon = min(abs(dLonAll));
        if minDLon < maxRes
            minDLon = maxRes;
        end
        if minDLat < maxRes
            minDLat = maxRes;
        end

        %define points to use and adjust ending bounds:
        stepR = min([minDLat,minDLon]);
        disp(ln80(['A stepsize of ' num2str(stepR) ' will be used for resampling the ensemble grids.']));
        nLonR = round(abs(lonAll(end)-lonAll(1))/stepR)+1;
        lonAll(end) = lonAll(1) + (nLonR - 1)*stepR;
        nLatR = round(abs(latAll(1)-latAll(end))/stepR)+1;
        latAll(1) = latAll(end) + (nLatR - 1)*stepR;

        %Define common grid:
        lonR = linspace(lonAll(1),lonAll(end),nLonR);
        latR = linspace(latAll(end),latAll(1),nLatR);
        latR = fliplr(latR)';

        %Interpolate all data to common grid:
        for ii = 1 : nFiles
            data{kk,ii} = area_int_2D(lon{kk,ii},lat{kk,ii},data{kk,ii},lonR,latR);
            hdr{kk,ii} = ESRI_hdr(lonR,latR, 'cor');
        end
    end

    szData = size(data{kk,1});
    %Calculate average and standard deviation at each point:
    dataAvg = reshape(nanmean(reshape([data{kk,:}],numel(data{kk,1}),[]),2),szData(1),szData(2));
    dataStd = reshape(nanstd(reshape([data{kk,:}],numel(data{kk,1}),[]),2),szData(1),szData(2));

    %%Write ensemble mean and standard deviation to file:
    %Define output path:
    dirOut{kk} = fullfile(path{kk,nFiles},'ensemble');

    %Define names:
    [~,nm,~] = fileparts(fileUse{kk,1});
    indRoot = regexpi(nm,'\d{1}');
    if ~isempty(indRoot)
        nm = nm(1:indRoot(end));
    end
    nm = [nm, '_ensemble'];
    fileAvg = fullfile(dirOut{kk}, [nm, '_avg.asc']);
    fileStd = fullfile(dirOut{kk}, [nm, '_std.asc']);
    disp(ln80(['The output files for the ensemble mean and standard ' ...
        'deviation are labeled ' char(39) fileAvg char(39) ' and ' ...
        char(39) fileStd char(39) ', respectively.' char(10)]));

    %Write files:
    write_ESRI_v4(dataAvg, hdr{kk,1}, fileAvg, 1);
    write_ESRI_v4(dataStd, hdr{kk,1}, fileStd, 1);
end
