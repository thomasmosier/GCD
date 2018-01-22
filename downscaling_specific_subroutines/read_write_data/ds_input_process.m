function [fldsDs, fldsBc] = ds_input_process(fldsTVar, varDs, dsPeriod)

%Find indice of field to downscale:
indFldDs = nan(0,1);
indFldBc = nan(0,1);


for ii = 1: numel(fldsTVar(:))
    if regexpbl(fldsTVar{ii}, {'sim', varDs}, 'and')
       indFldDs(end+1) = ii;
    elseif regexpbl(fldsTVar{ii}, 'ref')
        indFldBc(end+1) = ii;
    end
end
% if numel(fldsTVar(:)) > 0 && numel(fldsTVar(:)) < 3
%     for ii = 1: numel(fldsTVar(:))
%         if regexpbl(fldsTVar{ii}, 'sim')
%            indFldDs = ii;
%         elseif regexpbl(fldsTVar{ii}, 'ref')
%             indFldBc = ii;
%         end
%     end
% else
%     error('methodDirect:multTimeInv','Only one time invariant field expected.');
% end

%Extract fields for indicies
if all(~isnan(indFldDs))
    fldsDs = fldsTVar(indFldDs);
else
    fldsDs = cell.empty;
end

if all(~isnan(indFldBc))
    fldsBc = fldsTVar(indFldBc);
else
    fldsBc = cell.empty;
end

%Convert to cell (if character)
if ischar(fldsDs)
    fldsDs = {fldsDs};
end

if ischar(fldsBc)
    fldsBc = {fldsBc};
end

%Determine which field should actually be downscaled
if numel(fldsDs(:)) > 1
    indChck = min([numel(dsPeriod), 4]);
    for ii = numel(fldsDs(:)) : -1 : 1
        if ~regexpbl(fldsDs{ii}, dsPeriod(indChck))
            fldsDs(ii) = [];
        end
    end
end

if numel(fldsDs(:)) > 1
    warning('dsInputProcessing:multDsFlds',[num2str(numel(fldsDs(:))) ...
        ' fields were identified for downscaling. This is not expected and may cause an error.']);
elseif numel(fldsDs(:)) == 0
    error('dsInputProcessing:noDsFlds','No fields to downscale were found.');
end