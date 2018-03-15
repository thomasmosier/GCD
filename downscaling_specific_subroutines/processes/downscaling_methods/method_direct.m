function method_direct(sPath, sDs)


%Load time invariant data:
sTInv = ds_ld_fields(sPath, sDs.fldsTInv, sDs.lonDs, sDs.latDs, nan(1,2), nan(1,1));
if numel(sTInv(:)) == 1
    latOut = sTInv{1}.(sDs.varLat)(sTInv{1}.(sDs.varLat) >= min(sDs.latDs) & sTInv{1}.(sDs.varLat) <= max(sDs.latDs));
    lonOut = sTInv{1}.(sDs.varLon)(sTInv{1}.(sDs.varLon) >= min(sDs.lonDs) & sTInv{1}.(sDs.varLon) <= max(sDs.lonDs));
    
    [lonMeshOut, latMeshOut] = meshgrid(lonOut, latOut);
else
   error('methodDirect:multTimeInv','Only one time invariant field expected.');
end

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
    [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, 'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample);
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
    
    %Create input grid
    [lonMeshIn, latMeshIn] = meshgrid(sTVar{sDs.indDs}.(sDs.varLon), sTVar{sDs.indDs}.(sDs.varLat));
    
    %Initialize output structure:
    sTsOut = sTVar{sDs.indDs};
        sTsOut.lat = latOut;
        sTsOut.lon = lonOut;
        sTsOut.time = sTVar{sDs.indDs}.time;
        sTsOut.(sDs.varDs) = nan([numel(sTVar{sDs.indDs}.time), size(lonMeshOut)], 'single');

    %%ITERATE OVER timesteps:
    for jj = 1 : numel(sTsOut.date(:,1))
        if regexpi(sDs.intrp,'pchip')
            sTsOut.(sDs.varDs)(jj,:,:) = ...
                PCHIP_2D(lonMeshIn, latMeshIn, squeeze(sTVar{sDs.indDs}.(sDs.varDs)(jj,:,:)), lonMeshOut, latMeshOut);
        else
            sTsOut.(sDs.varDs)(jj,:,:) = ...
                interp2(lonMeshIn, latMeshIn, squeeze(sTVar{sDs.indDs}.(sDs.varDs)(jj,:,:)), lonMeshOut, latMeshOut, sDs.intrp);
        end
    end
    clear jj
    
    %Write inputs and outputs:
    if ~isempty(sDs.wrtOut)
%         ds_wrt_inputs(sTVar, sDs, sPath);
        
        if mnthCurr < 10
            strMnth = ['0' num2str(mnthCurr)];
        else
        	strMnth = num2str(mnthCurr);
        end
        strData = sPath.(['nm_' sDs.fldsTVar{indDsIn}]){1};
        
        for kk = 1 : numel(sDs.wrtOut(:))
            if regexpbl(sDs.wrtOut{kk}, {'ds','ts'}, 'and') || regexpbl(sDs.wrtOut{kk}, {'downscale','ts'}, 'and')
                %Make file and path names:
                fileCurr = [sDs.varDs '_' sDs.timestep '_' strData '-ds-' sDs.region '_intrp-' sDs.intrp '_' ...
                    num2str(min(sDs.yrsDs)) strMnth '01' '-' ...
                    num2str(max(sDs.yrsDs)) strMnth num2str(eomday(max(sDs.yrsDs), mnthCurr))];

                ds_wrt_outputs(sTsOut, 'intrp_ts', sDs, sPath, 'file', fileCurr, 'yrs', sDs.yrsDs);
            elseif regexpbl(sDs.wrtOut{kk}, {'ds','clim'}, 'and')
                %Create interpolated climatology
                sDsClm = sTsOut;
                    sDsClm.(sDs.varDs) = nanmean(sTsOut.(sDs.varDs)(sTsOut.date(:,1) >= min(sDs.yrsDs) & sTsOut.date(:,1) <= max(sDs.yrsDs),:,:), 1);
                    if isfield(sDsClm, 'time')
                        sDsClm.time = nan;
                    end
                    if isfield(sDsClm, 'date')
                        sDsClm.('date') = nan(1,2);
                    end
 
                %Make file and path names:
                fileCurr = [sDs.varDs '_' strData '-intrpclim-' sDs.region '_intrp-' sDs.intrp '_' ...
                    num2str(min(sDs.yrsDs)) 'thru' ...
                    num2str(max(sDs.yrsDs)) '_' strMnth];
                ds_wrt_outputs(sDsClm, 'intrp_clim', sDs, sPath, 'folder', 'intrp_clim', 'file', fileCurr); 
            end
        end
        clear kk
    end
    
    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes.' char(10)]);
end