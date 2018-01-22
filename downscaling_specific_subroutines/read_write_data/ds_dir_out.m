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

function foldOut = ds_dir_out(sPath, indSimCurr, foldSubRt)


%Find Name of time series directory:
fieldsPath = fieldnames(sPath);
nNms = numel(fieldsPath);
indTs = nan(nNms, 1);
indProj = nan(nNms, 1);
indSim = nan(nNms, 1);
for ii = 1 : nNms
    if regexpbl(fieldsPath{ii}, 'ts')
        indTs(ii) = ii;
    end
    if regexpbl(fieldsPath{ii}, 'proj')
        indProj(ii) = ii;
    end
    if regexpbl(fieldsPath{ii}, 'sim')
        indSim(ii) = ii;
    end
end
indTs(isnan(indTs)) = [];
indProj(isnan(indProj)) = [];
indSim(isnan(indSim)) = [];

if ~isempty(indProj)
    indUse = intersect(indTs,indProj);
elseif ~isempty(indTs)
    indUse = indTs;
else
   error('dsDirOut:noTs','No time-series folder was found. Therefore no output directory has been created.');
end

if ~isempty(indUse)
    if numel(indUse) > 1
        indUseTemp = intersect(indUse, indSim);
        if ~isempty(indUseTemp)
           indUse = indUseTemp(1); 
        end
    end
    indUse = indUse(1);
else
    error('dsDirOut:noTsUse','No useable time-series folder was found. Therefore no output directory has been created.');
end


%Find/Define output folder
foldRtCurr = sPath.(fieldsPath{indUse});

if iscell(foldRtCurr) 
    if numel(foldRtCurr) > 1
        if indSimCurr <= numel(foldRtCurr)
           foldRt = foldRtCurr{indSimCurr};
        else
            foldRt = foldRtCurr{1};
            warning('dsDirOut:tsInd',['The root output directory is being set to ' ...
                foldRt ' because the time-series simulation indice seems to be incorrect.'])
        end
    else
        foldRt = char(foldRtCurr);
    end
else
    foldRt = char(foldRtCurr);
end

[foldOutRt, fileLd, ~] = fileparts(foldRt); 

foldOut = fullfile(foldOutRt, foldSubRt);

%Check for simulation name (add if found):
indSimUnd = regexpi(fileLd,'_');
if ~isempty(indSimUnd) && length(indSimUnd) > 1
    nmSim = fileLd(indSimUnd(2)+1:indSimUnd(end)-1);
    
    if numel(nmSim) > 12
       nmSim = ''; 
    end
else
    nmSim = '';
end
if ~isempty(nmSim)
    foldOut = [foldOut '_' nmSim];
end


%Ensure output directory unique (doesn't write over existing data):
cntrDir = 0;
while exist(foldOut, 'dir')
    cntrDir = cntrDir + 1;
    if cntrDir == 1
        foldOut = [foldOut,'_' num2str(cntrDir)];
    else
        indC = regexpi(foldOut,'_');
        foldOut = [foldOut(1:indC(end)), num2str(cntrDir)];
    end
end

%Make directory:
if ~exist(foldOut, 'dir')
   mkdir(foldOut); 
end

%Display name of output path:
disp(['The output path is ' char(39) foldOut char(39) '.' char(10)]);
