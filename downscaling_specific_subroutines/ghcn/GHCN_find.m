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

function useData = GHCN_find(crd, dates, metVar, dataCon)



%This script loads all lat, lon, and station ID in from from GHCN *.inv 
%file. It then selects stations in specified bounding box using their lat
%and lon and returns the corresponding station ID.  It then loads the GHCN
%*.dat file, picks out the desired stations using their ID and desired 
%range of years, and returns the corresponding time-series data within a 
%cell array.

%coordinates = [lonW, lonE, latS, latN]
%dates = [starting year, ending year]
%metVar = 'pre' or 'tmn'
%dataCon = 'adj' or 'non'

%useData is a cell-array of size n by 5, where n is the number of stations
%within the bounding box defined by 'coordinates.'  {i,1} is the GHCN
%station ID, {i,2} is the station latitude, {i,3} is the station longitude,
%{i,4} is the station elevation, and {i,5} is itself a cell-array that 
%includes the time-series data.  {i,5} has dimensions m by 13, where m is
%the number of years of data present for the station (within dates 
%specified).  {j,1} is the year and {j,2:13} are the monthly values for the
%given year, given in mm in the case of precipitation and deg Celcius for 
%temperature.
lonW = crd(1);
lonE = crd(2);
latS = crd(3);
latN = crd(4);

startYear = dates(1);
endYear = dates(2);

if regexpi(metVar,'pre')    %Presumes all GHCN data are in subfolder 'GHCN data' within script folder.
    fidGhcnStn = fopen(which('v2.prcp.inv'));   %Stn ID (i11), Stn Name (a20), country (a10), latitude (f7.2), longitude (f8.2), and elevation (m) (i5)
    if ~isempty(regexpi(dataCon, 'adj'))
        fidGhcnData = fopen(which('v2.prcp_adj'));
    elseif regexpbl(dataCon,{'org','orig','non'})
        fidGhcnData = fopen(which('v2.prcp'));
    else
        error('Unknown data control type.');
    end
elseif regexpi(metVar,'tmp')
    if regexpi(dataCon, 'adj')
        fidGhcnData = fopen(which('ghcnm.tavg.v3.2.0.20120913.qca.dat'));
        fidGhcnStn = fopen(which('ghcnm.tavg.v3.2.0.20120913.qca.inv'));%Stn ID, Lat, Lon, Elevation
    elseif regexpbl(dataCon,{'org','orig','non'})
        fidGhcnData = fopen(which('ghcnm.tavg.v3.2.0.20120913.qcu.dat'));
        fidGhcnStn = fopen(which('ghcnm.tavg.v3.2.0.20120913.qcu.inv'));%Stn ID, Lat, Lon, Elevation
    else
        error('Unknown data control type.');
    end
else
   error('Invalid meteorological variable chosen.'); 
end

%PROCESS ALL STATION METADATA:
stnCellRaw = textscan(fidGhcnStn,'%s','Delimiter','');
fclose(fidGhcnStn);
stnMatTemp = cell2mat(stnCellRaw{1});
    clear stnCellRaw
    stnMat = nan(length(stnMatTemp),4); %Initialize
stnMat(:,1) = str2num(stnMatTemp(:,1:11));
    
if regexpi(metVar,'pre')      
    stnMat(:,2) = str2num(stnMatTemp(:,43:49));
    stnMat(:,3) = str2num(stnMatTemp(:,50:57));
    stnMat(:,4) = str2num(stnMatTemp(:,59:62));
elseif regexpi(metVar,'tmp') 
    stnMat(:,2) = str2num(stnMatTemp(:,13:20));
    stnMat(:,3) = str2num(stnMatTemp(:,22:30));
    stnMat(:,4) = str2num(stnMatTemp(:,32:37));
end
    clear stnMatTemp

%PROCESS ALL STATION DATA:   
dataCellRaw = textscan(fidGhcnData,'%s','Delimiter','');
fclose(fidGhcnData);
dataMatTemp = cell2mat(dataCellRaw{1});
    clear dataCellRaw
    dataMat = nan(length(dataMatTemp),14);   %Initialize
dataMat(:,1) = str2num(dataMatTemp(:,1:11));

if regexpi(metVar,'pre')   
    dataMat(:,2) = str2num(dataMatTemp(:,13:16));
    jj = 17;
    dj = 5;
    dataScale = 10;  %Data is given in tenths of milimeters
elseif regexpi(metVar,'tmp') 
    dataMat(:,2) = str2num(dataMatTemp(:,12:15));
    jj = 20;
    dj = 8;
    dataScale = 100;    %Temperature is given in hundreths of degrees celcius
end
   
for ii = 1 : 12
    dataMat(:,2+ii) = str2num(dataMatTemp( : , jj : jj+4 ))/dataScale;
    jj = jj + dj;
end  
    clear dataMatTemp
dataMat(dataMat == -9999/dataScale) = nan;

if regexpi(metVar,'pre')  
dataMat(dataMat == -8888/dataScale) = 0;
end

indRow = find( stnMat(:,2) > latS & stnMat(:,2) < latN & stnMat(:,3) > lonW & stnMat(:,3) < lonE);

useStn = stnMat(indRow,1:4);

useData = cell(length(useStn),5);   %Initialize

kk = 1; %Counter used in following loop
for ii = 1 : length(useStn(:,1))
    if ~isempty(dataMat(dataMat(:,1) == useStn(ii,1) & dataMat(:,2) >= startYear & dataMat(:,2) <= endYear , 2:end) )   %Ensures station has data for years of interest.
        useData{kk,1} = useStn(ii,1); %Station ID
        useData{kk,2} = useStn(ii,2);   %Latitude
        useData{kk,3} = useStn(ii,3);   %Longitude
        useData{kk,4} = useStn(ii,4);   %Elevation
        useData{kk,5} = dataMat(dataMat(:,1) == useStn(ii,1) & dataMat(:,2) >= startYear & dataMat(:,2) <= endYear,2:end); %Year, 12 months of data
        kk = kk + 1;
    end
end
useData(kk:end,:) = [];




