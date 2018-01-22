function [fldsTVar, fldsTInv] = ds_input_time_var(sPath)


%Find time variant and invariant fields:
fldsLd = fieldnames(sPath);

for ii = numel(fldsLd) : -1 : 1
    if strcmpi(fldsLd{ii}(1:2), 'nm')
        fldsLd(ii) = [];
    end
end
clear ii

fldsTInv = fldsLd;
fldsTVar = fldsLd;
for ii = numel(fldsLd) : -1 : 1
    if regexpbl(fldsLd{ii}, 'output')
        fldsTInv(ii) = [];
        fldsTVar(ii) = [];
    elseif regexpbl(fldsTInv{ii}, {'ts', 'clim', 'cm'}) && ~regexpbl(fldsTInv{ii}, 'pts')
        fldsTInv(ii) = [];
    else
        fldsTVar(ii) = [];
        
%         if regexpbl(fldsTInv{ii}, 'grid')
%            indUnd = regexpi(fldsTInv{ii}, '_');
%            fldsTInv{ii} = [fldsTInv{ii}(1:indUnd(end)) varElev];
%         end
    end
end
clear ii