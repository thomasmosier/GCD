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



%Warning: if a current file has the same name as the desired output, the
%existing file will be overwritten.

%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules

%%INNERWORKINGS (DO NOT MODIFY)
%Find directy of ESRI gridded ASCII files to convert:
%User Select folder path:
uiwait(msgbox(sprintf(['Select the folder of climate date in ESRI '...
    'ASCII format to convert to NetCDF format. \n']), ...
    '(Click OK to Proceed)','modal'));
dirData = uigetdir(['Select the directory containing '... 
    'the climate data in ESRI ASCII format.']);

%Find root filename and date range:
fileTemp = dir([dirData filesep '*.asc']);
filesTS = cell(length(fileTemp),1); 
for ii = 1 : length(fileTemp)
    filesTS(ii) = cellstr(fileTemp(ii).name);
end

%%Identify time-series elements present:
aryTime = NaN(length(filesTS),3);   %Initialize array to use

%Put all time-series elements into array:
for ii = 1 : length(filesTS)
    indSep = regexpi(filesTS{ii},'_|\.');
    aryTime(ii,1:2) = [str2double(filesTS{ii}(indSep(end-2)+1:indSep(end-1)-1)), ... 
        str2double(filesTS{ii}(indSep(end-1)+1:indSep(end)-1))];
end

%Sort array of time-series elements:
%Ensure Year comes first:
if length(aryTime(1,1)) < length(aryTime(1,2)) 
	aryTime = [aryTime(:,2), aryTime(:,1)];
end

%Sort by year, then month:
[aryTime, iTime] = sortrows(aryTime,[1,2]);
filesTS = filesTS(iTime);

%Determine number of days in each month within arytime:
aryTime(:,3) = eomday(aryTime(:,1),aryTime(:,2));

%%Write NetCDF File:
%Identify data type:
if ~isempty(regexpi(filesTS{1},'pre'))
    clmVarLng = 'precipitation';
    clmVar = 'pre';
elseif ~isempty(regexpi(filesTS{1},'tmx')) || ~isempty(regexpi(filesTS{1},'tasmax'))
    clmVarLng = 'maximum temperature';
    clmVar = 'tmx';
elseif ~isempty(regexpi(filesTS{1},'tmn')) || ~isempty(regexpi(filesTS{1},'tasmin'))
    clmVarLng = 'minimum temperature';
    clmVar = 'tmn';
elseif ~isempty(regexpi(filesTS{1},'tmp')) || ~isempty(regexpi(filesTS{1},'tas')) %'tas' must come after 'tasmax' and 'tasmin' to eliminate possibility of a false positive
    clmVarLng = 'mean temperature';  
    clmVar = 'tmp';
end

%Load example gridded data:
[dataEx, hdrESRI] = read_ESRI(fullfile(dirData,filesTS{1}));

%Initialize file:
pathNC = fullfile(dirData,[filesTS{1}(1:indSep(end-2)-1), '_', ...
    num2str(aryTime(1,1)), 'thru' num2str(aryTime(end,1)), '.nc']);
if exist(pathNC,'file') ~= 0
    ovrwrt = input(['A file with the name ' char(39) pathNC char(39) ...
        ' already exists.' char(10) ...
        'Would you like to overwrite this file? (Y or N) \n'],'s');
    if ~isempty(regexpi(ovrwrt,'y'))
        disp(['The file ' pathNC ' is being overwritten.'])
        delete(pathNC);
    else
        disp(['You chose to not overwrite the file ' pathNC ...
            '.  The script will now quit.']);
        quit;
    end
end

nccreate(pathNC,'time','Dimensions',{'time' length(aryTime(:,1))}, ...
    'Format','netcdf4', 'Datatype', 'single');
nccreate(pathNC,'latitude', 'Dimensions',{'lat' hdrESRI(2)}, ...
    'Format','netcdf4', 'Datatype', 'double');
nccreate(pathNC,'longitude','Dimensions',{'lon' hdrESRI(1)}, ...
    'Format','netcdf4', 'Datatype', 'double');
if ~isempty(regexpi(clmVar,'pre'))
    dataTyp = 'uint16';
else
    dataTyp = 'single';
end
nccreate(pathNC, clmVar, 'Dimensions', ...
    {'time' length(aryTime(:,1)) 'lat' hdrESRI(2) 'lon' hdrESRI(1)}, ...
    'Format','netcdf4', 'Datatype', dataTyp);

%Define attributes to be written
metaTime = {...   
    'bounds', 'time_bnds'; ...
    'units', 'days since 1900-01-01'; ...
    'calendar', 'gregorian'; ...
    'axis', 'T'; ...
    'long_name', 'time'; ...
    'standard_name', 'time'; ...
    '_CoordinateAxisType', 'Time'};
metaLon = {...
    'bounds', 'lon_bnds'; ...
    'units', 'degrees_east'; ...
    'axis', 'X'; ...
    'long_name', 'longitude'; ...
    'standard_name', 'longitude'; ...
    '_CoordinateAxisType', 'Lon'};
metaLat = {...
    'bounds', 'lat_bnds'; ...     
    'units', 'degrees_north'; ...
    'axis', 'Y'; ...
    'long_name', 'latitude'; ...
    'standard_name', 'latitude'; ...
    '_CoordinateAxisType', 'Lat'};

metaData = {...
        'standard_name', ''; ...
        'long_name', ''; ...
        'comment', ''; ...
        'units', ''; ...
        'cell_methods', 'time interval: days to middle of month (rounded); includes leap day)'; ...
        'cell_measures', 'area: geographic coordinates; coordinate is middle of each cell'; ...
        'history', ['converted ESRI ASCII formated data in ' dirData ' on ' date ';  Example filename: ' filesTS{1}]; ...
        'missing_value', num2str(hdrESRI(6))};

if ~isempty(regexpi(clmVar,'pre'))
    metaData{1,2} = 'precipitation_total'; 
    metaData{2,2} = 'Precipitation';
    metaData{3,2} = 'includes both liquid and solid phases from all types of clouds';
    metaData{4,2} = 'mm';
elseif ~isempty(regexpi(clmVar,'tmp')) || ~isempty(regexpi(clmVar,'tmn')) || ~isempty(regexpi(clmVar,'tmx'))
    if max(max(dataEx)) > 100
        metaData{4,2} = 'Kelvin';
    else
        metaData{4,2} = 'Celsius';
    end
    
    if ~isempty(regexpi(clmVar,'tmp'))
        metaData{1,2} = 'air_temperature'; 
        metaData{2,2} = 'Near-Surface Air Temperature';
        metaData{3,2} = 'monthly mean of the daily-mean near-surface air temperature.';
    elseif ~isempty(regexpi(clmVar,'tmn'))
        metaData{1,2} = 'air_temperature'; 
        metaData{2,2} = 'Daily Minimum Near-Surface Air Temperature';
        metaData{3,2} = 'monthly mean of the daily-minimum near-surface air temperature.';
    elseif ~isempty(regexpi(clmVar,'tmx'))
        metaData{1,2} = 'air_temperature'; 
        metaData{2,2} = 'Daily Maximum Near-Surface Air Temperature';
        metaData{3,2} = 'monthly mean of the daily-maximum near-surface air temperature.';
    end
end

strSrc = input(['Enter any source information to write into the ' ...
    'netCDF attributes list (e.g. low-res source dataset and version, '...
    'method notes, etc.):' char(10)],'s');

attGen = {...
    'title', ''; ...
    'institution', 'Data held at Oregon State University.'; ...
    'source', strSrc; ...
    'history', ['Converted on ' date ' using scripts created by Thomas Mosier (mosier.thomas@gmail.com) for use with associated downscaling package.']; ...
    'references', 'Information on the data is available at http://globalclimatedata.org'; ...
    'comment', 'Data restrictions: for academic research use only. Contact Thomas Mosier for details'; ...
    'reference', 'Mosier, T.M., Hill, D.F., Sharp, K.V., 2013. 30-Arcsecond monthly climate surfaces with global land coverage. International Journal of Climatology. DOI: 10.1002/joc.3829'; ...
    'contact', 'Thomas Mosier <mosier.thomas@gmail.com>'};

currYrTemp = date;
currYr = str2double(currYrTemp(end-3:end));
if aryTime(end,1) < currYr && ( (isempty(regexpi(filesTS{1},'gcm')) && isempty(regexpi(filesTS{1},'rcp'))) || ~isempty(regexpi(filesTS{1},'delta')) )
   attGen{1,2} = 'Delta downscaled monthly hindcast data';
else 
   attGen{1,2} = 'Downscaled GCM data';
end

%Initialize waitbar:
warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
hWait = waitbar(0,['The ESRI ASCII data in ' dirData ...
    ' is currently being processed.  '...
    'The first ASCII file will be converted to NetCDF format shortly.'],...
    'CreateCancelBtn','button_callback');

%Write global attributes:
for ii = 1 : length(attGen(:,1))
    ncwriteatt(pathNC,'/',attGen{ii,1},attGen{ii,2}) %'/' indicates global attribute
end

%Write time variable attributes:
for ii = 1 : length(metaTime(:,1))
    ncwriteatt(pathNC,'time',metaTime{ii,1},metaTime{ii,2}) %'/' indicates global attribute
end

%Write longitude variable attributes:
for ii = 1 : length(metaLon(:,1))
    ncwriteatt(pathNC,'longitude',metaLon{ii,1},metaLon{ii,2}) %'/' indicates global attribute
end

%Write latitude variable attributes:
for ii = 1 : length(metaLat(:,1))
    ncwriteatt(pathNC,'latitude',metaLat{ii,1},metaLat{ii,2}) %'/' indicates global attribute
end

%Write data variable attributes:
for ii = 1 : length(metaData(:,1))
    ncwriteatt(pathNC,clmVar,metaData{ii,1},metaData{ii,2}) %'/' indicates global attribute
end

%%Create and Write Variables to NetCDF File:
%Time:
vecTime = days_between([1900,1,1], [aryTime(:,1), aryTime(:,2), round(aryTime(:,3)/2)], 0);
ncwrite(pathNC,'time',vecTime);

%Longitude:
vecLon = hdrESRI(3) + 0.5*hdrESRI(5) : hdrESRI(5) : hdrESRI(3) + (hdrESRI(1)-0.5)*hdrESRI(5);
vecLon = vecLon';
ncwrite(pathNC,'longitude',vecLon);

%Latitude:
vecLat = hdrESRI(4) + (hdrESRI(2)-0.5)*hdrESRI(5) : -hdrESRI(5) : hdrESRI(4) + 0.5*hdrESRI(5);
vecLat = vecLat';
ncwrite(pathNC,'latitude',vecLat);

%Data:
for ii = 1 : length(filesTS(:))
    dataCurr = NaN([1, hdrESRI(2), hdrESRI(1)]);
    [dataCurr(1,:,:), ~] = read_ESRI(fullfile(dirData,filesTS{ii}));
    ptStrt = [ii,1,1];
%     stride = NaN(length(vecTime(:)),length(vecLat(:)),length(vecLon(:)));
%     stride(ii,:,:) = 1; 
%     ncwrite(pathNC,clmVar,dataCurr,ptStrt,stride);
    ncwrite(pathNC,clmVar,dataCurr,ptStrt);
    
    %Update waitbar:
    waitbar(ii/length(filesTS(:)), hWait, ['The file ' ...
        char(filesTS{ii}) ' has been written to ' char(pathNC) '.']);
end

delete(hWait);
warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
disp(['The script has finished.  ' num2str(length(filesTS(:))) ...
    ' ESRI formatted ASCII data files were converted to netCDF ' ...
    'format and placed in ' pathNC '.']);

% %For Testing and Diagnostics: 
% ncdisp(pathNC)
% for ii = 1 : length(filesTS(:))
%     dataSampl = ncread(pathNC,clmVar,[ii 1 1],[1 hdrESRI(2) hdrESRI(1)]);
%     dataSampl = squeeze(dataSampl);
%     [dataAct, ~] = read_ESRI(fullfile(dirData,filesTS{ii}));
%     if ~isequal(dataSampl,dataAct)
%         disp(['The data in ' fullfile(dirData,filesTS{ii}) ' is not equal to that in the ' num2str(ii) 'time element of ' pathNC]);
%     else
%         disp(['The ' num2str(ii) ' elements are equal'])
%     end
% end