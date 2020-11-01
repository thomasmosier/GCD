function sOut = stn_avg_bc(sStnRec, sSim, sSim2Bc, varDs)


varDate = 'date';

if regexpbl(varDs, 'pr')
    type = 'mult';
elseif regexpbl(varDs, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
    type = 'add';
end

%Find time-series elements to compare during bias correction:
%For historic model:
[indSimH, ~] = ismember(sSim.(varDate)(:,1:end), sStnRec.(varDate)(:,1:end), 'rows');
[indRefH, ~] = ismember(sStnRec.(varDate)(:,1:end), sSim.(varDate)(:,1:end),  'rows');

indSimH = find(indSimH ~= 0);
indRefH = find(indRefH ~= 0);

if numel(indSimH) ~= numel(indSimH)
    error('stnAvgBc:diffNumber',['Two different numbers of time-series ' ...
        'elements were selected. This is unexpected.']);
end

sOut = sSim;
if numel(indSimH) == 0 || numel(indRefH) == 0 
    warning('stnAvgBc:noStnRec',['No station records were found for times '...
        'overlapping with the simulation. Therefore the output will not be edited.']);
    return
else
    %Calculate bias between gridded reference (ERA) and station
    %observations (This is applied later to ERA and GCM before writing
    if strcmpi(type, 'mult')
        stnBias = nanmean(squeeze(sSim.(varDs)(indSimH,:,:))) / nanmean(sStnRec.(varDs));

        sOut.(varDs) = sSim2Bc.(varDs) / stnBias;
    elseif strcmpi(type, 'add')
        stnBias = nanmean(squeeze(sSim.(varDs)(indSimH,:,:))) - nanmean(sStnRec.(varDs));

        sOut.(varDs) = sSim2Bc.(varDs) - stnBias;
    else
        error('stnAvgBc:unknownBcType',[type ' is an unknown bias correction method.']);
    end
end