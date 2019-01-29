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

function ds_engine_v10(sPath, sMeta)


%Initialize downscaling struct to store lots of info that is standardized
%between methods
sDs = struct;

sDs.varDs = sMeta.varOut{sMeta.indDs};


[sDs.fldsTVar, sDs.fldsTInv] = ds_input_time_var(sPath);

sDs.varLat = 'latitude';
sDs.varLon = 'longitude';

sDs.latDs = sMeta.latOut{sMeta.indRegion};
sDs.lonDs = sMeta.lonOut{sMeta.indRegion};

if isfield(sMeta, 'indSim')
    sDs.indSim = sMeta.indSim;
end

if isfield(sMeta, 'nSim')
    sDs.nSim = sMeta.nSim;
end
    
    
    
if isfield(sMeta, 'indRegion')
    sDs.region = sMeta.region{sMeta.indRegion};
end

if isfield(sMeta, 'wrtOut')
    sDs.wrtOut = sMeta.wrtOut;
elseif isfield(sMeta, 'wrtData')
    sDs.wrtOut = sMeta.wrtData;
else
    warning('methodDirect:noOutput',['No field indicating the fields '...
        'to write to file was found. Therefore, only the downscaled output will be written.']);
    sDs.wrtOut = {'downscale_ts'};
end

if isfield(sMeta, 'wrtPrec')
    sDs.wrtPrec = sMeta.wrtPrec;
elseif isfield(sMeta, 'prec')
    sDs.wrtPrec = sMeta.prec;
else
    %Define precision of output files:
    if regexpbl(sDs.varDs, 'pre') || strcmpi(sDs.varDs, 'pr')
        sDs.wrtPrec = 1;   %Write precipitation as integers
    elseif regexpbl(sDs.varDs, {'tav', 'tmn', 'tmean', 'tmp', 'tas', 'tmx', 'tmax', 'tmn', 'tmin'})
        sDs.wrtPrec = 2;   %Write temperature with one decimal  
    else
        sDs.wrtPrec = 2;
        warning('dsEngine:unknownVar',[sDs.varDs ' was selected as the '...
            'climate variable to processed, which has not been programmed for. This may cause the program to crash.']);
    end

%     warning('methodDirect:noWrtPrec',['No output precision was specified. '...
%         'It is being set to ' num2str(sDs.wrtPrec) '.']);
end


if isfield(sMeta, 'period')
    sDs.period = sMeta.period;
else
    error('dsEngine:unknownPeriod',['No field was found denoting the ' ...
        'downscaling period (either projection or historical).']);
end


if isfield(sMeta, 'timestep')
    sDs.timestep = sMeta.timestep;
elseif isfield(sMeta, 'timeStep')
    sDs.timestep = sMeta.timeStep;
else
    sDs.timestep = 'unknown';
    warning('dsEngine:noTimestep',['No timestep field found. Therefore, '...
        'the timestep is being set to ' sDs.timestep]);
end


if isfield(sMeta, 'method')
    sDs.method = sMeta.method;
else
    error('dsEngine:unknownMethod',['No field was found denoting the ' ...
        'downscaling method was found.']);
end

if isfield(sMeta, 'intrp')
    sDs.intrp = sMeta.intrp;
elseif isfield(sMeta, 'interp')
    sDs.intrp = sMeta.interp;
else
    error('dsEngine:unknownIntrpMethod',['No field was found denoting the ' ...
        'downscaling interpolation method was found.']);
end

if isfield(sMeta, 'upscale')
   sDs.resample = sMeta.upscale;
elseif isfield(sMeta, 'resample')
   sDs.resample = sMeta.resample;
else
    error('dsEngine:unknownUpscaleMethod',['No field was found denoting the ' ...
        'downscaling resampling/upscaling method.']);
end

if isfield(sMeta, 'uiInput')
    sDs.methodspec = sMeta.('uiInput'){sMeta.indDs};
elseif isfield(sMeta, 'uiinput')
    sDs.methodspec = sMeta.('uiinput'){sMeta.indDs};
end
    
    


if isfield(sMeta, 'wrtFormat') 
    sDs.wrtFmt = sMeta.wrtFormat;
elseif isfield(sMeta, 'wrtTyp')
    sDs.wrtFmt = sMeta.wrtTyp;
else
    sDs.wrtFmt = 'asc';
    warning('methodDirect:noWrtFormat',['No output file format was specified. '...
        'Therefore, it is being set to ' sDs.wrtFmt '.']);
end

%Override default units:
if isfield(sMeta,'unitsUse')
    sDs.units = sMeta.unitsUse;
elseif isfield(sMeta,'unitsuse')
    sDs.units = sMeta.unitsuse;
elseif isfield(sMeta,'units')
    sDs.units = sMeta.units;
else
    sDs.units = cell(0,2);
end

%Months out:
if isfield(sMeta, 'mnthsOut')
    sDs.mnthsDs = sMeta.mnthsOut;
elseif isfield(sMeta, 'mnths')
    sDs.mnthsDs = sMeta.mnths;
elseif isfield(sMeta, 'mnthsLd')
    sDs.mnthsDs = sMeta.mnthsLd;
else
    sDs.mnthsDs = (1:12);
    warning('methodDirect:noMnthsFound',['No field indicating the months '...
        'to process was found. Therefore, all months will be used.']);
end

sDs.yrsLd = [min([min(sMeta.yrsOut), min(sMeta.yrsBase)]) max([max(sMeta.yrsOut), max(sMeta.yrsBase)])];

sDs.yrsBase = [min(sMeta.yrsBase), max(sMeta.yrsBase)];
sDs.yrsDs = [min(sMeta.yrsOut), max(sMeta.yrsOut)];


if isfield(sMeta, 'yrsBc')
    sDs.yrsBc = sMeta.yrsBc;
    sDs.yrsLd = [min([min(sDs.yrsLd), min(sDs.yrsBc)]) max([max(sDs.yrsLd), max(sDs.yrsBc)])];
else
    sDs.yrsBc = sDs.yrsDs;
end

%Find indice of field to downscale:
[sDs.fldsDs, sDs.fldsBc] = ds_input_process(sDs.fldsTVar, sDs.varDs, sDs.period);

if numel(sDs.fldsDs(:)) == 0 
    error('dsEngine:noDsFlds','No downscaing fields were found. One should be found.')
elseif numel(sDs.fldsDs(:)) > 1
    if regexpbl(sDs.method, '2var')
        for ii = numel(sDs.fldsDs) : -1 : 1
            if regexpbl(sDs.fldsDs{ii}, '4bc')
                sDs.fldsDs(ii) = [];
            end
        end
        
        if numel(sDs.fldsDs(:)) > 1
            error('dsEngine:multDsFlds2Var','Multiple downscaing fields were found. One should be found.');
        end
    else
        error('dsEngine:multDsFlds','Multiple downscaing fields were found. One should be found.');
    end
end

%Set indice of data to downscale
indDs = find(strcmpi(sDs.fldsTVar, sDs.fldsDs{1}) == 1);
if ~isempty(indDs)
    sDs.indDs = indDs;
else
   sDs.indDs = nan; 
end

%Set indice of data to use for bias correction (if present)
if isfield(sDs, 'fldsBc') && ~isempty(sDs.fldsBc(:))
    sDs.indRef = nan(numel(sDs.fldsBc(:)),1);
    for ii = 1 : numel(sDs.fldsBc(:))
        sDs.indRef(ii) = find(strcmpi(sDs.fldsTVar, sDs.fldsBc{ii}) == 1);
    end
    clear ii
else
    sDs.indRef = nan;
end

sDs.pathOut = fullfile(sPath.output, 'output');
if ~exist(sDs.pathOut, 'dir')
    mkdir(sDs.pathOut);
end

if isfield(sMeta, 'shepard') && sMeta.shepard == 1
	disp(['NOTE: Bias Correction using inverse distance weighting ' ...
        'has been selected which adds about a factor of 4 to the ' ...
        'processing time.']);
end
if regexpbl(sMeta.wrtTyp, 'asc')
	disp(['NOTE: A significant amount of time '...
        'is spent writing data to ASCII files.  The speed can be '...
        'reduced substantially by limiting the number of outputs.']);
end
disp(blanks(1));   %Blank line


%Check that GPCC is only used for precipitation:
fldsPath = fieldnames(sPath);
for ii = 1 : numel(sPath)
    if regexpbl(sPath.(fldsPath{ii}),'GPCC') && ~regexpbl(sDs.varDs, 'pre')
        error('dsEngine:GPCC_tmn',['GPCC is only valid as a low-res' ...
            ' precipitation source.  They do not produce temperature data.']);
    end
end
clear ii

%Create sPath version in which only current simulation is present
sPathCurr = sPath;
fldsPath = fieldnames(sPath);
for ii = numel(fldsPath) : -1 : 1
    if regexpbl(fldsPath{ii}, {'nm'})
        fldsPath(ii) = []; 
    end
end

for ii = 1 : numel(fldsPath)
    if iscell(sPathCurr.(fldsPath{ii})) && numel(sPathCurr.(fldsPath{ii})(:)) == sDs.nSim
        sPathCurr.(fldsPath{ii}) = sPathCurr.(fldsPath{ii}){sDs.indSim};
        fldNm = ['nm_' fldsPath{ii}];
        if isfield(sPathCurr, fldNm) && iscell(sPathCurr.(fldNm)) && numel(sPathCurr.(fldNm)(:)) == sDs.nSim
            sPathCurr.(fldNm) = sPathCurr.(fldNm){sDs.indSim};
        end
    end
end
if iscell(sPathCurr.('output')) && numel(sPathCurr.('output')(:)) == sDs.nSim
    sPathCurr.('output') = sPathCurr.('output'){sDs.indSim};
end


%%IMPLEMENT DOWNSCALING METHOD
disp(['NOTE: This program does not freeze, but it can take several ' ...
    'minutes for certain aspects of the script to complete,' char(10) ' depending ' ...
    'on the method, geographic domain, and time period.']);
if regexpbl(sMeta.method, 'delta')
    method_delta(sPathCurr, sDs);
% elseif regexpbl(sMeta.method, 'bcsd')
%     method_bcsd(sPathCurr, sDs)
elseif regexpbl(sMeta.method, 'direct')
    method_direct(sPathCurr, sDs);
elseif regexpbl(sMeta.method, 'extract')
    method_extract(sPathCurr, sDs);
elseif regexpbl(sMeta.method, {'bc','only'})
    method_bc_only(sPathCurr, sDs);
elseif regexpbl(sMeta.method, 'tlapse')
    if ~regexpbl(sDs.varDs, {'tmp','tmx','tmn','tas','tave','tavg','tmin','tmax'})
        error('dsEngine:wrongVar',['method_tlapse is only designed to work with temperature variables (current variable is ' sDs.varDs ')']);
    end
    
    method_tlapse(sPathCurr, sDs);
elseif regexpbl(sMeta.method, 'pw') 
    if ~regexpbl(sDs.varDs, 'pre') && ~strcmpi(sDs.varDs, 'pr')
        error('dsEngine:wrongVar',['method_pw is only designed to work with precipitation variables (current variable is ' sDs.varDs ')']);
    end
    method_pw(sPathCurr, sDs);
elseif regexpbl(sMeta.method, 'roe') 
    if ~regexpbl(sDs.varDs, 'pre') && ~strcmpi(sDs.varDs, 'pr')
        error('dsEngine:wrongVar',['method_pw is only designed to work with precipitation variables (current variable is ' sDs.varDs ')']);
    end
    method_roe(sPathCurr, sDs);
else
    error('dsEngine:unknownMethod', [char(39) sMeta.method char(39) ' has not been programmed for.']);
end


