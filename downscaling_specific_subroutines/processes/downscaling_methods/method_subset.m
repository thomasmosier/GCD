function method_subset(sPath, sDs)


%Load time invariant data:
% sTInv = ds_ld_fields(sPath, sDs.fldsTInv, sDs.lonDs, sDs.latDs, nan(1,2), nan(1,1));
% if numel(sTInv(:)) == 1
%     latOut = find( sTInv{1}.(sDs.varLat) >= min(sDs.latDs) & sTInv{1}.(sDs.varLat) <= max(sDs.latDs));
%     lonOut = find( sTInv{1}.(sDs.varLon) >= min(sDs.lonDs) & sTInv{1}.(sDs.varLon) <= max(sDs.lonDs));
%     
%     [lonMeshOut, latMeshOut] = meshgrid(lonOut, latOut);
% else
%    error('methodDirect:multTimeInv','Only one time invariant field expected.');
% end

indDsIn = sDs.indDs;

%Loop over months (or lump all months together if daily)
if regexpbl(sDs.timestep, {'month', 'mnth'})
    sDs.nLp = numel(sDs.mnthsDs);
elseif regexpbl(sDs.timestep, {'day', 'daily'})
    sDs.nLp = 1;
else
    warning('methodDirect:unknownTimestep',['The timestep is ' ...
        sDs.timestep ', which has not been programmed for.']);
end


tic;
for ii = 1 : sDs.nLp
    if regexpbl(sDs.timestep, {'month', 'mnth'})
        mnthCurr = sDs.mnthsDs(ii);
        mnthDisp = mnth_str(ii);
    elseif regexpbl(sDs.timestep, {'day', 'daily'})
        mnthCurr = nan;
        mnthDisp = 'annual';
    end
    
    %%LOAD INPUTS:
        indDsOrg = sDs.indDs;
    [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, ...
        'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample, 'units', sDs.units);
    %Update 'indDsIn' based on change in indDs
        sDs.indDs = indNew(1);
        sDs.indRef = indNew(2:end);
        indDelta = (sDs.indDs - indDsOrg);
        indDsIn = indDsIn + indDelta;
        
    %BIAS CORRECT INPUT
    if ~isempty(sDs.indRef) && ~isnan(sDs.indRef) 
        [sTVar{end+1}, bcMethod] = bc_switch(sTVar, indDsIn, sDs);
        sDs.fldsTVar{numel(sTVar(:))} = [sDs.fldsTVar{indDsIn} '_bc-' bcMethod];
        sDs.indDs = numel(sDs.fldsTVar);
%         sTVar{indDsIn} = bc_switch(sTVar, sDownscale, sMeta);
    else 
        sDs.indDs = indDsIn;
    end

    %Write inputs and outputs:
    if ~isempty(sDs.wrtOut)
%         ds_wrt_inputs(sTVar, sDs, sPath)
        ds_wrt_outputs(sTVar{sDs.indDs}, 'output', sDs, sPath)
    end
    
    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes.' char(10)]);
end