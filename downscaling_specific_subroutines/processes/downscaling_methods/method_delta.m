function method_delta(sPath, sDs)


varLon = 'longitude';
varLat = 'latitude';
varDate = 'date';

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

%Find indice of high-resolution reference climatology
sDs.indHrClm = nan;
for ii = numel(sDs.fldsTVar) : -1 : 1
    if regexpbl(sDs.fldsTVar{ii}, {'hr', 'clim', sDs.varDs}, 'and')
        sDs.indHrClm = ii;
    end
end

indDsIn = sDs.indDs;

%Set name of data product
strData = sPath.(['nm_' sDs.fldsTVar{indDsIn}]){1};
if numel(strData) > 10
   strData = 'output'; 
end

%Loop over months (or lump all months together if daily)
if regexpbl(sDs.timestep, {'month', 'mnth'})
    sDs.nLp = numel(sDs.mnthsDs);
else
    error('methodDelta:wrongTimestep',['The timestep is ' sDs.timestep ...
        ' but the Delta method should only be used with a monthly timestep.']);
end



%Determine delta type based on variable:
if regexpbl(sDs.varDs, 'pr')
    deltaType = 'mult';
elseif regexpbl(sDs.varDs, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
    deltaType = 'add';
else
    error('methodDelta:unknownVar',['Variable ' sDs.varDs ' has not '...
        'been programmed for. Determine which delta method to use.']);
end

tic;
for ii = 1 : sDs.nLp
    if regexpbl(sDs.timestep, {'month', 'mnth'})
        mnthCurr = sDs.mnthsDs(ii);
    end

    %%LOAD INPUTS:
        indDsOrg = sDs.indDs;
    [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, 'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample, 'frame', 2);
    %Update 'indDsIn' based on change in indDs
        sDs.indDs = indNew(1);
        sDs.indRef = indNew(2:end);
        indDelta = (sDs.indDs - indDsOrg);
        indDsIn = indDsIn + indDelta;

    if ~regexpbl(sTVar{indDsIn}.timestep, {'month', 'mnth'})
        error('methodDelta:dataWrongTimestep',['The data timestep is ' sTVar{indDsIn}.timestep ...
        ' but the Delta method should only be used with a monthly timestep.'])
    end

    %Get output grid:
    indLatOut = find( sTVar{sDs.indHrClm}.(sDs.varLat) >= min(sDs.latDs) & sTVar{sDs.indHrClm}.(sDs.varLat) <= max(sDs.latDs));
    indLonOut = find( sTVar{sDs.indHrClm}.(sDs.varLon) >= min(sDs.lonDs) & sTVar{sDs.indHrClm}.(sDs.varLon) <= max(sDs.lonDs));

    %Crop high-resolution reference climatology:
    sTVar{sDs.indHrClm}.(sDs.varLat) = sTVar{sDs.indHrClm}.(sDs.varLat)(indLatOut);
    sTVar{sDs.indHrClm}.(sDs.varLon) = sTVar{sDs.indHrClm}.(sDs.varLon)(indLonOut);
    
    if numel(size(sTVar{sDs.indHrClm}.(sDs.varDs))) == 3
        sTVar{sDs.indHrClm}.(sDs.varDs) = sTVar{sDs.indHrClm}.(sDs.varDs)(:,indLatOut,indLonOut);
    elseif numel(size(sTVar{sDs.indHrClm}.(sDs.varDs))) == 2
        sTVar{sDs.indHrClm}.(sDs.varDs) = sTVar{sDs.indHrClm}.(sDs.varDs)(indLatOut,indLonOut);
    else
       error('methodDelta:unknownHrClimSize',['The high-resolution '...
           'reference climatology array has ' ...
           num2str(numel(size(sTVar{sDs.indHrClm}.(sDs.varDs)))) ...
           ' dimensions. 2 or 3 are expected.']) 
    end
    
    [lonMeshOut, latMeshOut] = meshgrid(sTVar{sDs.indHrClm}.(sDs.varLon), sTVar{sDs.indHrClm}.(sDs.varLat));
    
    %BIAS CORRECT INPUT
    if isfield(sDs, 'indRef') && all(~isnan(sDs.indRef))
        for ll = 1 : numel(sDs.indRef)
            if ~regexpbl(sTVar{sDs.indRef(ll)}.timestep, {'month', 'mnth'})
                error('methodDelta:bcWrongTimestep',['The data timestep is ' sTVar{sDs.indRef(ll)}.timestep ...
                    ' but the Delta method should only be used with a monthly timestep.'])
            end
        end
        clear ll
    
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
        sTsOut.(varLat) = sTVar{sDs.indHrClm}.(sDs.varLat);
        sTsOut.(varLon) = sTVar{sDs.indHrClm}.(sDs.varLon);
        sTsOut.time = sTVar{sDs.indDs}.time;
        sTsOut.(sDs.varDs) = nan([numel(sTVar{sDs.indDs}.time), size(lonMeshOut)], 'single');

    %Create simulation climatology
    indSimClm = find(sTVar{sDs.indDs}.(varDate)(:,1) >= min(sDs.yrsBase) & sTVar{sDs.indDs}.(varDate)(:,1) <= max(sDs.yrsBase));
    nYrsBase = (max(sDs.yrsBase) - min(sDs.yrsBase) + 1);
    if numel(indSimClm) < nYrsBase
        error('methodDelta:missingSimClimYrs',[num2str(numel(indSimClm)) ...
           ' time indices were found to produce the simulation climatology (' ...
           num2str(sTVar{sDs.indDs}.(varDate)(indSimClm(1),1)) '-' ...
           num2str(sTVar{sDs.indDs}.(varDate)(indSimClm(end),1)) '), but ' ...
           num2str(nYrsBase) ' are needed (' num2str(min(sDs.yrsBase)) '-' ...
           num2str(max(sDs.yrsBase)) '). One cause of this is if there is if data in input folder has a temporal discontinuity.']);
    end
    sSimClm = sTVar{sDs.indDs};
        sSimClm.(sDs.varDs) = squeeze(nanmean(sTVar{sDs.indDs}.(sDs.varDs)(indSimClm,:,:), 1));
        sSimClm.time = nan;
        sSimClm.date = [nan, mnthCurr, nan];

    %Initialize input anomaly:
    sAnomIn = sTVar{sDs.indDs};
        sAnomIn.(sDs.varDs) = nan(size(sTVar{sDs.indDs}.(sDs.varDs)), 'single');
        
    %Initialize output anomaly:
    sAnomOut = sTsOut;
        sAnomOut.(sDs.varDs) = nan(size(sTsOut.(sDs.varDs)), 'single');
    
    %%ITERATE OVER timesteps:
    for jj = 1 : numel(sTsOut.time)
        %Calculate anomaly:
        if regexpbl(deltaType, 'mult')
            sAnomIn.(sDs.varDs)(jj,:,:) = squeeze(sTVar{sDs.indDs}.(sDs.varDs)(jj,:,:)) ./ sSimClm.(sDs.varDs);
        elseif regexpbl(deltaType, 'add')
            sAnomIn.(sDs.varDs)(jj,:,:) = squeeze(sTVar{sDs.indDs}.(sDs.varDs)(jj,:,:))  - sSimClm.(sDs.varDs);
        end
        
        %Interpolate anomaly:
        if regexpi(sDs.intrp,'pchip')
                warning('off', 'PCHIP_2D:NaN');
            sAnomOut.(sDs.varDs)(jj,:,:) = ...
                PCHIP_2D(lonMeshIn, latMeshIn, squeeze(sAnomIn.(sDs.varDs)(jj,:,:)), lonMeshOut, latMeshOut);
                warning('on', 'PCHIP_2D:NaN');
        else
            sAnomOut.(sDs.varDs)(jj,:,:) = ...
                interp2(lonMeshIn, latMeshIn, squeeze(sAnomIn.(sDs.varDs)(jj,:,:)), lonMeshOut, latMeshOut, sDs.intrp);
        end
        
        %Scale anomaly by high resolution climatology:
        if regexpbl(deltaType, 'mult')
            sTsOut.(sDs.varDs)(jj,:,:) = squeeze(sAnomOut.(sDs.varDs)(jj,:,:)) .* squeeze(sTVar{sDs.indHrClm}.(sDs.varDs));
        elseif regexpbl(deltaType, 'add')
            sTsOut.(sDs.varDs)(jj,:,:) = squeeze(sAnomOut.(sDs.varDs)(jj,:,:))  + squeeze(sTVar{sDs.indHrClm}.(sDs.varDs));
        end
    end
    clear jj
    
    %Write inputs and outputs:
    if ~isempty(sDs.wrtOut)
        ds_wrt_inputs(sTVar, sDs, sPath);
        
        if mnthCurr < 10
            strMnth = ['0' num2str(mnthCurr)];
        else
        	strMnth = num2str(mnthCurr);
        end
        
        %Write outputs
        for kk = 1 : numel(sDs.wrtOut(:))
            if regexpbl(sDs.wrtOut{kk}, {'ds','ts'}, 'and') || regexpbl(sDs.wrtOut{kk}, {'downscale','ts'}, 'and')
                %Make file and path names:
                fileCurr = [sDs.varDs '_' sDs.timestep '_' strData '-ds-' sDs.region '_delta_intrp-' sDs.intrp '_' ...
                    num2str(min(sDs.yrsDs)) strMnth '01' '-' ...
                    num2str(max(sDs.yrsDs)) strMnth num2str(eomday(max(sDs.yrsDs), mnthCurr))];

                ds_wrt_outputs(sTsOut, 'delta', sDs, sPath, 'file', fileCurr, 'yrs', sDs.yrsDs);
            elseif regexpbl(sDs.wrtOut{kk}, {'hr','anom'}, 'and')
                ds_wrt_outputs(sAnomOut, 'hranom', sDs, sPath, 'folder', fullfile('hranom'), 'yrs', sDs.yrsDs);
            elseif regexpbl(sDs.wrtOut{kk}, {'lr','anom'}, 'and')
                ds_wrt_outputs(sAnomIn, 'lranom', sDs, sPath, 'folder', fullfile('lranom'), 'yrs', sDs.yrsDs);
            elseif regexpbl(sDs.wrtOut{kk}, {'in','clim'}, 'and')
                ds_wrt_outputs(sSimClm, 'inputclim', sDs, sPath, 'folder', fullfile('inputclim'));
            elseif regexpbl(sDs.wrtOut{kk}, {'ref','clim'}, 'and')
                ds_wrt_outputs(sSimClm, 'refclim', sDs, sPath, 'folder', fullfile('refclim'));    
            elseif regexpbl(sDs.wrtOut{kk}, {'ds','clim'}, 'and')
                sDsClm = sSimClm;
                    sDsClm.(varLon) = sTsOut.(varLon);
                    sDsClm.(varLat) = sTsOut.(varLat);
                    if isfield(sDsClm, 'time')
                        sDsClm.time = nan;
                    end
                    if isfield(sDsClm, varDate)
                        sDsClm.(varDate) = nan(1,2);
                    end
        
                sDsClm.(sDs.varDs) = nanmean(sTsOut.(sDs.varDs)(sTsOut.date(:,1) >= min(sDs.yrsDs) & sTsOut.date(:,1) <= max(sDs.yrsDs),:,:), 1);
                ds_wrt_outputs(sDsClm, 'dsclim', sDs, sPath, 'folder', fullfile('dsclim')); 
            end
        end
    end
    
    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnth_str(ii) ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes.' char(10)]);
end