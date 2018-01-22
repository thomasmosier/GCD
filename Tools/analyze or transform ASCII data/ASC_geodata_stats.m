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



%Script finds info from ASCII files in selected folder and tabulated info
%based upon defined queries:




%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules

%Get folder with ts data:
uiwait(msgbox(sprintf(['Select the directory containing '... 
    'ESRI gridded ASCII monthly data. \n']),...
    '(Click OK to Proceed)','modal'));
foldTs = uigetdir(['Select the directory containing '... 
    'ESRI gridded ASCII monthly data.']);
disp([foldTs ' has been selected as the folder of monthly ESRI ASCII files.']);

%Find files in folder:
filesTemp = dir(fullfile(foldTs,'*.asc'));
gridFiles = struct2cell(filesTemp);
gridFiles = gridFiles(1,:);

%Get file with reference data:
uiwait(msgbox(sprintf(['Select a reference file in the same format as '...
    'the time-series data (same geographical extent) \n']),...
    '(Click OK to Proceed)','modal'));
[fileRef, pathRef, ~] = uigetfile({'*.asc';'*.txt'},['Select a reference file in the same format as '...
    'the time-series data (same geographical extent)']);
disp([fullfile(pathRef, fileRef) ' has been selected as the reference file.']);

%Initiate waitbar:
hWait = waitbar(0,'Data processing has begun.');

%Load reference file:
[dataRef,hdrRef,metaRef] = read_ESRI(fullfile(pathRef, fileRef));

%Find mean of reference data:
refMn = mean(reshape(dataRef,1,[]));
refStd = std(reshape(dataRef,1,[]));

dataRefGTmn = zeros(size(dataRef));
dataRefStdMid = zeros(size(dataRef));
dataRefGTmn(dataRef >= refMn) = 1;
dataRefStdMid(dataRef >= refMn - refStd & dataRef <= refMn + refStd) = 1;

%Write these to files:
write_ESRI_v4(dataRefGTmn, hdrRef, fullfile(pathRef, [fileRef(1:end-4) 'GTmean.asc']), 0);
write_ESRI_v4(dataRefStdMid, hdrRef, fullfile(pathRef, [fileRef(1:end-4) 'midSTD.asc']), 0);

%%Create Output directory:
dirOutput = mk_output_dir(fullfile(foldTs, gridFiles{1}), 'geodata_stats');
disp(['The output files will be written to: ' dirOutput]);
indUnd = regexpi(gridFiles{1},'_');
fileOut = fullfile(dirOutput, gridFiles{1}(1:indUnd(end-1)-1));

%Create comma seperated file with stats:
cellHdr = {'Year','Month','Avg (entire)','Avg (>=Mean alt)', 'Avg (<Mean alt)', 'Avg (Middle Std)'};
strHdr = blanks(0);

for ii = 1 : numel(cellHdr)
%     strHdr = [strHdr char(cellHdr{ii}) ', '];
    if ii ~= numel(cellHdr)
        strHdr = [strHdr char(cellHdr{ii}) ', '];
    else
        strHdr = [strHdr char(cellHdr{ii}) char(10)];
    end
end

%Write Header:
fileID = fopen(fileOut,'w');
fwrite(fileID, strHdr, 'char');

%Initialize output matrix:
nFiles = numel(gridFiles);
matWrt = nan(nFiles,numel(cellHdr));

for ii = 1 : nFiles
    %Update waitbar:
    warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    fracComplt = ii/nFiles;
    waitbar(fracComplt, hWait, ...
        ['The file ' char(39) char(strrep(gridFiles{ii},'_','\_')) ...
        char(39) ' is being read and processed.']);
    warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    
    %Find time corresponding to file:
    indFUnd = regexpi(gridFiles{ii},'_');
    indExt = regexpi(gridFiles{ii},'\.');
    matWrt(ii,1) = str2double(gridFiles{ii}(indFUnd(end-1)+1:indFUnd(end)-1));
    matWrt(ii,2) = str2double(gridFiles{ii}(indFUnd(end)+1:indExt(end)-1));
    [dataCurr,hdr,meta] = read_ESRI(fullfile(foldTs, gridFiles{ii}));
    [lat, lon] = ESRI_hdr2geo(hdr, meta);
    area = area_geodata(lon,lat);
    
    matWrt(ii,3) = nansum(nansum(dataCurr.*area)) / nansum(nansum(area));
    matWrt(ii,4) = nansum(nansum(dataCurr(dataRefGTmn == 1).*area(dataRefGTmn == 1))) / nansum(nansum(area(dataRefGTmn == 1)));
    matWrt(ii,5) = nansum(nansum(dataCurr(dataRefGTmn == 0).*area(dataRefGTmn == 0))) / nansum(nansum(area(dataRefGTmn == 0)));
    matWrt(ii,6) = nansum(nansum(dataCurr(dataRefStdMid == 1).*area(dataRefStdMid == 1))) / nansum(nansum(area(dataRefStdMid == 1)));
end


%Sort:
matWrt = sortrows(matWrt,[1 2]);

%Write data:
dlmwrite(fileOut,matWrt,'-append');

fclose(fileID);
delete(hWait);