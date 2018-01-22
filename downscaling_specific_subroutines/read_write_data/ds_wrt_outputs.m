function ds_wrt_outputs(sTsOut, type, sDs, sPath, varargin)

foldRt = '';
fileRt = '';
yrs = nan(1,2);
var = sDs.varDs;
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
       if strcmpi(varargin{ii}, 'file')
           fileRt = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'fold') || strcmpi(varargin{ii}, 'folder')
           foldRt = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'yrs') || strcmpi(varargin{ii}, 'years')
           yrs = varargin{ii+1};
       elseif strcmpi(varargin{ii}, 'var')
           var = varargin{ii+1};
       end
    end
end

%Determine if main output:
if regexpbl(type, {'output', 'anom', 'ds','delta','extract','analysis'}) || regexpbl(type, {'in', 'clim'}, 'and')
    blMainOut = 1;
else
    blMainOut = 0;
end


%Set folder:
if ~isempty(foldRt)
    foldCurr = fullfile(sDs.pathOut, foldRt);
elseif blMainOut == 1
    foldCurr = sDs.pathOut;
else
    foldCurr = fullfile(sDs.pathOut, type);
end

%Set file:
if ~isempty(fileRt)
    fileCurr = fileRt;
elseif blMainOut
    nmOut = ds_file_nm(sPath, sDs, sDs.indDs);
    
    fileCurr = [var, '_', nmOut, '_', type, '_', sDs.region];
else
    fileCurr = [var, '_', type, '_', sDs.region];
end

%Crop years (if option selected)
if all(~isnan(yrs))
    if isfield(sTsOut, 'date')
        indUse = find(sTsOut.date(:,1) >= min(yrs) & sTsOut.date(:,1) <= max(yrs));
        
        sTsOut.(var) = sTsOut.(var)(indUse,:,:);
        sTsOut.('date') = sTsOut.('date')(indUse,:);
        if isfield(sTsOut, 'time')
        sTsOut.('time') = sTsOut.('time')(indUse,:);
        end
    else
        warning('dsWrtOutputs:noDateFld','The time is not being cropped because the date field does not exist.')
    end
end

%Crop lat/lon
varLon = 'longitude';
varLat = 'latitude'; 
%Find indices for cropping output to specificed downscaling area:
indLonOut = find(min(sDs.lonDs) <= sTsOut.(varLon) & sTsOut.(varLon) <= max(sDs.lonDs));
indLatOut = find(max(sDs.latDs) >= sTsOut.(varLat) & sTsOut.(varLat) >= min(sDs.latDs));

if ~isempty(indLonOut) && ~isempty(indLatOut) 
    %Crop outut data:
    sTsOut.(var) = sTsOut.(var)( : , min(indLatOut) : max(indLatOut) , min(indLonOut) : max(indLonOut));
    sTsOut.(varLat) = sTsOut.(varLat)(indLatOut);
    sTsOut.(varLon) = sTsOut.(varLon)(indLonOut);
end

%Write to file
write_geodata_v2(fullfile(foldCurr, fileCurr), sTsOut, sDs.wrtPrec, sDs.wrtFmt, 'var', var);