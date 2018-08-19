function ds_wrt_outputs(sTsOut, typ, sDs, sPath, varargin)

foldRt = '';
fileRt = '';
yrs = nan(1,2);
var = sDs.varDs;
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
        if ischar(varargin{ii})
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
end

%Determine if main output:
if regexpbl(typ, {'output', 'ds','delta','extract','analysis'})
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
    foldCurr = fullfile(sPath.output, 'extras', typ);
    if ~exist(foldCurr, 'dir')
        mkdir(foldCurr);
    end
end

%Set file:
if ~isempty(fileRt)
    fileCurr = fileRt;
elseif blMainOut
    nmOut = ds_file_nm(sPath, sDs, sDs.indDs);
    
    fileCurr = [var, '_', nmOut, '_', typ, '_', sDs.region];
else
    fileCurr = [var, '_', typ, '_', sDs.region];
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


if iscell(fileCurr)
    if numel(fileCurr(:)) == 1
        blWrt1 = 1;
        fileCurr = fileCurr{1};
    else
        blWrt1 = 0;
        if numel(fileCurr(:)) ~= numel(sTsOut.(var)(:,1,1))
            error('dsWrtOutputs:lengthMismatch',['There are ' ...
                num2str(numel(fileCurr(:))) ' output files and ' ...
                num2str(numel(sTsOut.(var)(:,1,1))) ' time slices. The two should be the same.'])
        end
    end
else
    blWrt1 = 1;
end

%Write to file
if blWrt1 == 1
    write_geodata_v2(fullfile(foldCurr, fileCurr), sTsOut, sDs.wrtPrec, sDs.wrtFmt, 'var', var);
else
    strDisp = '';
    for ii = 1 : numel(fileCurr(:))
        sTsOutCurr = sTsOut;
        sTsOutCurr.(var) =  squeeze(sTsOutCurr.(var)(ii,:,:));
        write_geodata_v2(fullfile(foldCurr, fileCurr{ii}), sTsOut, sDs.wrtPrec, sDs.wrtFmt, 'var', var, strDisp);
        
        if ii == 1
            strDisp = 'no_disp';
        end
    end
end