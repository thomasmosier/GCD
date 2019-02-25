function method_direct(sPath, sDs)

varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';


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
fldsIn = sDs.fldsTVar;

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
    [sTVar, ~, sDs.fldsTVar] = ds_ld_fields(sPath, fldsIn, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, ...
        'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample, 'frame', 2, 'units', sDs.units);
%     [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, 'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample);
    %Update indices of data to use
    for zz = 1 : numel(sDs.fldsTVar)
        if regexpbl(sDs.fldsTVar{zz}, 'sim')
            sDs.indDs = zz;
        elseif regexpbl(sDs.fldsTVar{zz}, 'ref')
            sDs.indRef = zz;
        end
    end

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

    
    %Define unique year-month combinations to loop over:
    datesYrMnth = unique(sTVar{sDs.indDs}.(varDate)(:,1:2), 'rows');
    %Keep only those specified for downscaling:
    datesYrMnth = datesYrMnth(datesYrMnth(:,1) >= min(sDs.yrsDs) & datesYrMnth(:,1) <= max(sDs.yrsDs), : );

    
    %Loop over output files:
    for ll = 1 : numel(datesYrMnth(:,1))   
        %Initialize grids for current month:
        indOutCurr = find(ismember(sTVar{sDs.indDs}.(varDate)(:,1:2), datesYrMnth(ll,1:2), 'rows'));
        
        %Initialize output structure:
        sTsOut = sTVar{sDs.indDs};
            sTsOut.(varLat) = latOut;
            sTsOut.(varLon) = lonOut;
            sTsOut.time = sTVar{sDs.indDs}.time(indOutCurr);
            sTsOut.(varDate) = round(sTVar{sDs.indDs}.(varDate)(indOutCurr,:));
            sTsOut.(sDs.varDs) = nan([numel(indOutCurr), size(lonMeshOut)], 'single');
        
        
        %%ITERATE OVER timesteps:
        for jj = 1 : numel(sTsOut.(varDate)(:,1))
            if regexpi(sDs.intrp,'pchip')
                sTsOut.(sDs.varDs)(jj,:,:) = ...
                    PCHIP_2D(lonMeshIn, latMeshIn, squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indOutCurr(jj),:,:)), lonMeshOut, latMeshOut);
            else
                sTsOut.(sDs.varDs)(jj,:,:) = ...
                    interp2(lonMeshIn, latMeshIn, squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indOutCurr(jj),:,:)), lonMeshOut, latMeshOut, sDs.intrp);
            end
        end
        clear jj

        
        %Write inputs and outputs:
        if ~isempty(sDs.wrtOut)
            nmStrCurr = strrep(['nm_' sDs.fldsTVar{indDsIn}], 'projhist', 'proj');
            strData = sPath.(nmStrCurr){1};

            for kk = 1 : numel(sDs.wrtOut(:))
                if regexpbl(sDs.wrtOut{kk}, {'ds','ts'}, 'and') || regexpbl(sDs.wrtOut{kk}, {'downscale','ts'}, 'and')
                     %Make file and path names:
                    fileDs = ds_output_name(sDs, strData, mnthCurr, 'ds');

                    ds_wrt_outputs(sTsOut, 'delta', sDs, sPath, 'file', fileDs, 'yrs', sDs.yrsDs);
                elseif regexpbl(sDs.wrtOut{kk}, {'in','clim'}, 'and')
                    ds_wrt_outputs(sSimClm, 'inputclim', sDs, sPath);
                elseif regexpbl(sDs.wrtOut{kk}, {'ref','clim'}, 'and')
                    ds_wrt_outputs(sSimClm, 'refclim', sDs, sPath, 'folder'); 
                elseif regexpbl(sDs.wrtOut{kk}, {'ds','clim'}, 'and')
                    %Create interpolated climatology
                    sDsClm = sTsOut;
                        sDsClm.(sDs.varDs) = nanmean(sTsOut.(sDs.varDs)(sTsOut.(varDate)(:,1) >= min(sDs.yrsDs) & sTsOut.(varDate)(:,1) <= max(sDs.yrsDs),:,:), 1);
                        if isfield(sDsClm, 'time')
                            sDsClm.time = nan;
                        end
                        if isfield(sDsClm, varDate)
                            sDsClm.(varDate) = nan(1,2);
                        end

                    %Make file and path names:
                    fileClim = ds_output_name(sDs, strData, mnth, 'clim');
                    ds_wrt_outputs(sDsClm, 'dsclim', sDs, sPath, 'folder', 'dsclim', 'file', fileClim); 
                end
            end
            clear kk
        end  
    end
    clear ll
    
    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes.' char(10)]); 
end
clear ii