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


%Currently does not work - problem with .grb files?
%-Maybe incomplete path used when calling GCM to cdfToolbox?

%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules

yrs = [1950, 2000];
mnths = (1:12);
crd = [70, 80.1, 33.3, 37.2];
region = 'NPakistan';
metVar = 'pre';

%Locate MACA NETCDF File:
uiwait(msgbox(sprintf(['Select '...
        'the GCM netCDF of GRIB file for ' ...
        metVar '.' char(10)]), '(Click OK to Proceed)','modal'));
[fileGcm, pathGcm, ~] = uigetfile({'*.nc'; '*.grb'}, ...
    ['GCM ' metVar ' Time-Series Selection']);

strucHdr = struct('metVar', metVar, 'yrsOut', yrs, 'mnths', mnths, ...
    'crd',crd, 'region',region);
clear('yrs','mnths','crd','region');

fullpathGcm = fullfile(pathGcm, fileGcm);
    
[gcmTs, units, gcmLat, gcmLon, yrsOutCurr] = GCM_load(fullpathGcm, strucHdr, metVar, []);

hdrGcm = gcm_hdr(gcmLat, gcmLon);

%Create output directory:
dirOutput = mk_output_dir(pathGcm, [strucHdr.region, '_ESRI_ASCII_' char(metVar)]);

hWait = waitbar(0,['The GCM ' fullpathGcm ' is currently being processed.'...
    'The first ASCII file will be written shortly.']);
for jj = 1 : length(gcmTs(:,1,1))
    currDate = [ceil(jj/12) + yrsOutCurr(1) - 1, rem(jj,12)];   %[Year, month]
    if currDate(2) == 0 %'Rem' has a range of 0->11, where 0 is actually the 12th month
       currDate(2) = 12; 
    end

    fileRtInd = regexpi(fileGcm, '_');
    filename = [char(fileGcm(1:fileRtInd(end)-1)) '_' ...
        num2str(currDate(1)) '_' num2str(currDate(2)) '.asc'];
    fullFilepath = fullfile(dirOutput, filename);

    dataCurr = squeeze(gcmTs(jj,:,:));

    if regexpi(units,'kg m-2 s-1') %Units are a flux (i.e. per time), so convert
        daysCurr = eomday(currDate(1), currDate(2));
        dataCurr = dataCurr * (daysCurr*24*60*60);
    end

    %Add functionality to clip data using reference ASCII file
    if exist('dataClip','var') && ~isempty(dataClip)
        [dataRef, hdrRef] = read_ESRI(dataClip);

        XvectRef = hdrRef(3) + 0.5*hdrRef(5) : hdrRef(5) : ...
            hdrRef(3) + (hdrRef(1) - 0.5 )*hdrRef(5) ;
        YvectRef = ( hdrRef(2) - 0.5 )*hdrRef(5) + hdrRef(4) : ...
            -hdrRef(5) : hdrRef(4) + 0.5*hdrRef(5) ;

        [XRef,YRef] = meshgrid( XvectRef, YvectRef);

        XvectGcm = hdrGcm(3) + 0.5*hdrGcm(5) : hdrGcm(5) : ...
            hdrGcm(3) + (hdrGcm(1) - 0.5 )*hdrGcm(5) ;
        YvectGcm = ( hdrGcm(2) - 0.5 )*hdrGcm(5) + hdrGcm(4) : ...
            -hdrGcm(5) : hdrGcm(4) + 0.5*hdrGcm(5) ;

        [XGcmMat,YGcmMat] = meshgrid( XvectGcm, YvectGcm);

        resampledRef = interp2(XRef,YRef,dataRef, XGcmMat,YGcmMat,'nearest'); 

        dataCurr( isnan(resampledRef) ) = NaN; 
    end

    write_ESRI_v3(dataCurr, hdrGcm, fullFilepath, 3);

    warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString');
    waitbar(jj/length(gcmTs(:,1,1)), hWait, ...
        ['The GCM file ' char(strrep(filename,'_','\_')) ' has been written.']);
    warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString');
end
delete(hWait);