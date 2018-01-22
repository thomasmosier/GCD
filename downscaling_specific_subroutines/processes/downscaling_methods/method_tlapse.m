function method_tlapse(sPath, sDs)

%Do not normalize grids:
strNorm = 'none';
% warning('methodTLapse:manualOpts','Set degN and elevHg in systematic/automatic way.')

strFigRes = '-r300';


%Make figures
blWrtFig = 1;

%Path for storing output figures:
if blWrtFig == 1
    foldFig = fullfile(sPath.output, 'figures');
    if ~exist(foldFig, 'dir')
        mkdir(foldFig)
    end
end
warning('off', 'MATLAB:LargeImage');

%SET VARIABLE NAMES (THESE ARE BASICALLY ARBITRARY):
varDate = 'date'; %Name to use for dates
varArea = 'area';
varLat = 'latitude';
varLon = 'longitude';

%General slopes/changes
varDeltaX = 'deltax'; %Name to use for x change
varDeltaY = 'deltay'; %Name to use for y change
varDeltaMag = 'deltamag'; %Name to use for combined gradient (sqrt[x^2+y^2])
varTLps = 'tlapse';

%Elevation
varDem = 'z'; %Data variable used in geopotential height field (also used to define elevation in other arrays)
varDemAnom = 'elevAnomMag'; %Variable to use for elevation anomaly
varDeltaYDemAnom = 'deltayElevAnom';
varDeltaXDemAnom = 'deltaxElevAnom';
varDeltaMagDemAnom = 'deltamagElevAnom';

if ~regexpbl(sDs.varDs, {'tas','tmp','tave','tmx','tmax','tmn','tmin'})
   error('methodTLapse:wrongVar',['The tlapse method is not designed to work with ' sDs.varDs]);
end


%Get method options:
if isfield(sDs, 'methodspec')
    lapseGrp = '';
    for ii = 1 : numel(sDs.methodspec(:,1))
        if regexpbl(sDs.methodspec{ii,1},{'lapse', 'type'})
            lapseGrp = sDs.methodspec{ii,2};
        end
    end
    clear ii
else
    error('methodTLapse:noMethodspec',['The field methodspec is required '...
        'to specificy method-specific options. This field was not found.']);
end

%Load time invariant data:
sTInv = ds_ld_fields(sPath, sDs.fldsTInv, sDs.lonDs, sDs.latDs, nan(1,2), nan(1,1));

%Find high-res output DEM:
indHrDem = nan;
indLrDem = nan;
if numel(sTInv(:)) > 0
    
    for ii = 1 : numel(sDs.fldsTInv)
       if regexpbl(sDs.fldsTInv{ii}, {'hr','dem'}, 'and') 
           indHrDem = ii;
       elseif regexpbl(sDs.fldsTInv{ii}, {'lr','dem'}, 'and')
           indLrDem = ii;
       end
    end
end
    
%Process high-res DEM
if ~isnan(indHrDem)
    sTInv{indHrDem}.(varDem) = single(squeeze(sTInv{indHrDem}.(varDem)));

    sTInv{indHrDem}.(varArea) = area_geodata(sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), 'c');
    
%     latOut = find( sTInv{indHrDem}.(sDs.varLat) >= min(sDs.latDs) & sTInv{indHrDem}.(sDs.varLat) <= max(sDs.latDs));
%     lonOut = find( sTInv{indHrDem}.(sDs.varLon) >= min(sDs.lonDs) & sTInv{indHrDem}.(sDs.varLon) <= max(sDs.lonDs));
%     
    [lonMeshHr, latMeshHr] = meshgrid(sTInv{indHrDem}.(sDs.varLon), sTInv{indHrDem}.(sDs.varLat));
else
   error('methodTlapse:noHrDem','No high-resolution DEM was found.');
end

%Process low-res DEM
if ~isnan(indLrDem)
    [lonMeshLr, latMeshLr] = meshgrid(sTInv{indLrDem}.(varLon), sTInv{indLrDem}.(varLat));
    sTInv{indLrDem}.(varDem) = squeeze(sTInv{indLrDem}.(varDem));
    sTInv{indLrDem}.(varArea) = area_geodata(sTInv{indLrDem}.(varLon), sTInv{indLrDem}.(varLat), 'c');

    %Calculate x,y changes in elevation
    [sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY)] = gradient_map(sTInv{indLrDem}.(varDem));
    [~, sTInv{indLrDem}.(varDeltaMag)] = cart2pol(sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY));
    
else
   error('methodTlapse:noLrDem','No low-resolution DEM was found.');
end

if ~isnan(indLrDem) && ~isnan(indHrDem)
    indDemIntrp = numel(sTInv(:)) + 1;
    sTInv{indDemIntrp} = sTInv{indHrDem};
    
    if regexpi(sDs.intrp,'pchip')
            warning('off', 'PCHIP_2D:NaN');
        sTInv{indDemIntrp}.(varDem) = ...
            PCHIP_2D(lonMeshLr, latMeshLr, squeeze(sTInv{indLrDem}.(varDem)), lonMeshHr, latMeshHr);
            warning('on', 'PCHIP_2D:NaN');
    else
        sTInv{indDemIntrp}.(varDem) = interp2(lonMeshLr, latMeshLr, ...
            squeeze(sTInv{indLrDem}.(varDem)), ...
            lonMeshHr, latMeshHr, sDs.intrp);
    end

    %Calculate difference between interpolated DEM and HR DEM:
    %This is used later on for determining lapse rate correct to apply
    sTInv{indDemIntrp}.(varDemAnom) = sTInv{indHrDem}.(varDem) - sTInv{indDemIntrp}.(varDem);

    %Calculate x,y changes in high-res elevation DEM anomaly (used for calculating alignment angle):
    [sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom)] ...
        = gradient_map(sTInv{indDemIntrp}.(varDemAnom));
    [~, sTInv{indDemIntrp}.(varDeltaMagDemAnom)] = cart2pol(...
        sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom));    
else
    error('methodTlapse:DemInputsNotFound','Either the low or high resolution DEM inputs have not been located.');
end


%Plot elevation histograms
if blWrtFig == 1
    hFig = figure('color', 'white','visible','off'); hold on; 
    hPClrY = pcolor(sTInv{indLrDem}.(varLon), sTInv{indLrDem}.(varLat), double(squeeze(sTInv{indLrDem}.(varDeltaY)))); colorbar; 
    hold off
    xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
    pathFig = fullfile(foldFig, 'lowres_deltaY');
    if exist([pathFig, '.fig'], 'file')
        delete([pathFig, '.fig'])
    end
    if exist([pathFig, '.png'], 'file')
        delete([pathFig, '.png'])
    end
    savefig(hFig, [pathFig '.fig']);
%         export_fig([pathPlot '.eps'],'-painters');
    export_fig([pathFig '.png'],'-painters',strFigRes);
    
    hFig = figure('color', 'white','visible','off'); hold on; 
    histogram(sTInv{indLrDem}.(varDeltaMag)(:));
    xlabel('Change Between Grid Cells (m)'); ylabel('Occurences');
    pathFig = fullfile(foldFig, 'low_res_elevation_change_mag');
    if exist([pathFig, '.fig'], 'file')
        delete([pathFig, '.fig'])
    end
    if exist([pathFig, '.png'], 'file')
        delete([pathFig, '.png'])
    end
    savefig(hFig, [pathFig '.fig']);
    export_fig([pathFig '.png'],'-painters',strFigRes);
    
        %Make figure:
    
    
%     %high-res elevation change magnitude
%     hFig = figure('color', 'white','visible','off'); hold on; 
%     histogram(sTInv{indHrDem}.(varDeltaMag)(:), 100);
%     xlabel('Change Between Grid Cells (m)'); ylabel('Occurences');
%     pathFig = fullfile(foldFig, 'High_res_elevation_change_mag');
%     if exist([pathFig, '.fig'], 'file')
%         delete([pathFig, '.fig'])
%     end
%     if exist([pathFig, '.png'], 'file')
%         delete([pathFig, '.png'])
%     end
%     savefig(hFig, [pathFig '.fig']);
%     export_fig([pathFig '.png'],'-painters',strFigRes);
    
    %Difference in SRTM and interpolated ERA DEMs
    hFig = figure('color', 'white','visible','off'); hold on; 
    histogram(sTInv{indDemIntrp}.(varDemAnom)(:), 100);
    xlabel('High-Res - Interpolated Low-Res DEMs (m)'); ylabel('Occurences');
    pathFig = fullfile(foldFig, 'High-Low_res_DEM_anom');
    if exist([pathFig, '.fig'], 'file')
        delete([pathFig, '.fig'])
    end
    if exist([pathFig, '.png'], 'file')
        delete([pathFig, '.png'])
    end
    savefig(hFig, [pathFig '.fig']);
    export_fig([pathFig '.png'],'-painters',strFigRes);
    
    %Magnitude of difference in SRTM and ERA DEMs
    hFig = figure('color', 'white','visible','off'); hold on; 
    histogram(sTInv{indDemIntrp}.(varDeltaMagDemAnom)(:), 100);
    xlabel('Mag. Changes in High and Low DEMs (m)'); ylabel('Occurences');
    pathFig = fullfile(foldFig, 'High-Low_res_DEM_anom_mag');
    if exist([pathFig, '.fig'], 'file')
        delete([pathFig, '.fig'])
    end
    if exist([pathFig, '.png'], 'file')
        delete([pathFig, '.png'])
    end
    savefig(hFig, [pathFig '.fig']);
%         export_fig([pathPlot '.eps'],'-painters');
    export_fig([pathFig '.png'],'-painters',strFigRes);
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
elseif regexpbl(sDs.timestep, {'day', 'daily'})
    sDs.nLp = 1;
else
    warning('methodDirect:unknownTimestep',['The timestep is ' ...
        sDs.timestep ', which has not been programmed for.']);
end

%Main processing loop:
tic;
for ii = 1 : sDs.nLp
    if regexpbl(sDs.timestep, {'month', 'mnth'})
        mnthCurr = sDs.mnthsDs(ii);
        mnthUse = mnthCurr;
        mnthDisp = mnth_str(ii);
    elseif regexpbl(sDs.timestep, {'day', 'daily'})
        mnthCurr = nan;
        mnthUse = (1:12);
        mnthDisp = 'annual';
    end
    
    if ~isequal(mnthUse, (1:12))
       error('methodTlapse:missingMonths', ['The method_tlapse function is '...
           'currently only designed to work when all 12 months are '...
           'simultaneously being downscaled.']); 
    end
    
    %%LOAD INPUTS:
        indDsOrg = sDs.indDs;
    [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, 'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample);
    %Update 'indDsIn' based on change in indDs
        sDs.indDs = indNew(1);
        sDs.indRef = indNew(2:end);
        indDelta = (sDs.indDs - indDsOrg);
        indDsIn = indDsIn + indDelta;
        
    if ~isequal(sTInv{indLrDem}.(varLon), sTVar{sDs.indDs}.(varLon))
        error('methodTLapse:londiff','The low-res DEM and simulation time-series longitude grids to not align.');
    elseif ~isequal(sTInv{indLrDem}.(varLat), sTVar{sDs.indDs}.(varLat))
        error('methodTLapse:latdiff','The low-res DEM and simulation time-series latitude grids to not align.');
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
    
%     %Create input grid
%     [lonMeshIn, latMeshIn] = meshgrid(sTVar{sDs.indDs}.(sDs.varLon), sTVar{sDs.indDs}.(sDs.varLat));
%     
%     %Initialize output structure:
%     sTsOut = sTVar{sDs.indDs};
%         sTsOut.lat = latOut;
%         sTsOut.lon = lonOut;
%         sTsOut.time = sTVar{sDs.indDs}.time;
%         sTsOut.(sDs.varDs) = nan([numel(sTVar{sDs.indDs}.time), size(lonMeshOut)], 'single');

    
    %Write inputs and outputs:
    if ~isempty(sDs.wrtOut) && ~isempty(sDs.indRef) && ~isnan(sDs.indRef) 
%         ds_wrt_inputs(sTVar, sDs, sPath)
        ds_wrt_outputs(sTVar{sDs.indDs}, 'bc', sDs, sPath)
    end
    
    
    %Create climatology(ies):
    sDs.indCm = numel(sTVar(:)) + 1;
    sTVar{sDs.indCm} = struct;
        sTVar{sDs.indCm}.(varLon) = sTVar{sDs.indDs}.(varLon);
        sTVar{sDs.indCm}.(varLat) = sTVar{sDs.indDs}.(varLat);
        sTVar{sDs.indCm}.(varDate) = [nan([numel(mnthUse),1]), mnthUse(:), nan([numel(mnthUse),1])];
        sTVar{sDs.indCm}.('time') = nan([numel(mnthUse),1]);
        sTVar{sDs.indCm}.(sDs.varDs) = nan([numel(mnthUse), size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:)))], 'single');
        sTVar{sDs.indCm}.(varDeltaX)   = sTVar{sDs.indCm}.(sDs.varDs);
        sTVar{sDs.indCm}.(varDeltaY)   = sTVar{sDs.indCm}.(sDs.varDs);
        sTVar{sDs.indCm}.(varDeltaMag) = sTVar{sDs.indCm}.(sDs.varDs);
        sTVar{sDs.indCm}.(varTLps)     = sTVar{sDs.indCm}.(sDs.varDs);
        
    %Calculate climatology:
    for mm = 1 : numel(mnthUse)
        indCurr = find(sTVar{sDs.indDs}.(varDate)(:,1) >= min(sDs.yrsBase) & sTVar{sDs.indDs}.(varDate)(:,1) <= max(sDs.yrsBase) & sTVar{sDs.indDs}.(varDate)(:,2) == mnthUse(mm));
        sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:) = squeeze(nanmean(sTVar{sDs.indDs}.(sDs.varDs)(indCurr,:,:), 1));
        
        if blWrtFig == 1
            hFig = figure('color', 'white','visible','off'); hold on; 
            hPClr = pcolor(sTVar{sDs.indCm}.(varLon), sTVar{sDs.indCm}.(varLat), double(squeeze(sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:)))); colorbar; 
            hold off
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [sDs.varDs '_clim_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
        %         export_fig([pathPlot '.eps'],'-painters');
            export_fig([pathFig '.png'],'-painters',strFigRes);
        end
    end
    clear mm
    
    
    %Calculate sptail changes in climatologies:
    for mm = 1 : numel(mnthUse)
        [sTVar{sDs.indCm}.(varDeltaX)(mm,:,:), sTVar{sDs.indCm}.(varDeltaY)(mm,:,:)] = gradient_map(squeeze(sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:)));

        qLr = alignment_factor(sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY), ...
            squeeze(sTVar{sDs.indCm}.(varDeltaX)(mm,:,:)), squeeze(sTVar{sDs.indCm}.(varDeltaY)(mm,:,:)));

        %Scale combined climatological change by alignment factor, q.
        [~, tempMagCmDs] = cart2pol(squeeze(sTVar{sDs.indCm}.(varDeltaX)(mm,:,:)), squeeze(sTVar{sDs.indCm}.(varDeltaY)(mm,:,:)));
        sTVar{sDs.indCm}.(varDeltaMag)(mm,:,:) = qLr.*tempMagCmDs;

        %Ratio between magnitude of gradient in climatology and topography is the lapse rate!
        sTVar{sDs.indCm}.(varTLps)(mm,:,:) = squeeze(sTVar{sDs.indCm}.(varDeltaMag)(mm,:,:))./sTInv{indLrDem}.(varDeltaMag);
    end
    
    %Calculate lapse rates:
    %Define lapse rate regions and calculate monthly lapse rates for each:
    if strcmpi(lapseGrp, '5NoSo')
        %Group areas as the following:
            %Northern portion of domain and steep
            %Southern portion of domain and steep
            %High elevation and flat
            %Northern portion of domain and flat
            %Southern portion of domain and flat
    
        warning('methodTLapse:lapserateFactors', ['The lapse rate parameters ' ...
            'degNo and elevHg are being set to values designed for the Himalaya ' ...
            'region. Consider changing for other regions.']);
        degNo = 35;
        elevHg = 2.5*10^3;  
        deltaSteep = 300;
%         degNo  = mean2d(sTInv{indLrDem}.(varLat));
%         elevHg = mean2d(sTInv{indLrDem}.(varDem));


        %ORIGINAL LAPSE RATES (CONSIDERS ASPECT; ALTERED VERSION DOES NOT):
        %Indices for all grid cells:
        indAll = (1:numel(sTInv{indLrDem}.(varDeltaMag)));
        %Grid cells that are steep:
        indSteep = find(sTInv{indLrDem}.(varDeltaMag) > deltaSteep);  
        %Grid cells that are flat:
        indFlat = setdiff(indAll, indSteep);
        %Grid cells in Northern portion of domain:
        indNo = find(latMeshLr(:) > degNo);  
        %Grid cells in Southern portion of domain:
        indSo = find(latMeshLr(:) <= degNo);
        %Grid cells that are in Southern portion of domain and South facing:
        indSoFc = intersect(indSo, find(sTInv{indLrDem}.(varDeltaY)(:) > 0));
        %Grid cells that are South facing and steep:
        indSoSteep = intersect(indSoFc, indSteep);
        %Grid cells that are steep, not South Facing, and in Northern part of
        %domain:
        indNoSteep = intersect(setdiff(indSteep, indSoSteep), indNo);
        %Grid cells that are high and flat:
        indHiFlat = intersect(find(sTInv{indLrDem}.(varDem)(:) > elevHg), indFlat); 
        %Grid cells that are in Northern portion of domain and flat/low elevation:
        indNoExt = setdiff(setdiff(indNo, indNoSteep), indHiFlat);
        %Grid cells that are in Southern portion of domain and flat/low elevation:
        indSoExt = setdiff(setdiff(indSo, indSoSteep), indHiFlat);
        %Indices to check:
        indChck = [indNoSteep; indSoSteep; indHiFlat; indNoExt; indSoExt];

                
        %Check to make sure that the five areas are unique and cover all indices:
        if numel(indChck) ~= numel(indAll) || numel(indChck) ~= numel(unique(indChck))
            error('ERA_ds:incorrectGrps',['The indices for the groups selected '...
                'are either non-unique or do not include all grid cells.']);
        end


        %Calculate monthly lapse rats for each group
        lrLpsMnth = nan([12, size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:)))], 'single');
        for zz = 1 : 12
            dataCurr = squeeze(sTVar{sDs.indCm}.(varTLps)(zz,:,:));
            
            %Calculate lapse rates
            lpRtNoSteep = nanmean(dataCurr(indNoSteep)); %Northern and steep
            lpRtSoSteep = nanmean(dataCurr(indSoSteep)); %Southern and steep
            lpRtHiFlat = nanmean(dataCurr(indHiFlat)); %High and flat
            lpRtNoFlat = nanmean(dataCurr(indNoExt)); %Extra Northern
            lpRtSoFlat = nanmean(dataCurr(indSoExt)); %Extra Southern
%             %Calculate lapse rates
%             lpRtNoSteep = nanmean(dataCurr(indNoSteep)); %Northern and steep
%             lpRtSoSteep = nanmean(dataCurr(indSoSteep)); %Southern and steep
%             lpRtHiFlat = nanmean(dataCurr(indHiFlat)); %High and flat
%             lpRtNoFlat = nanmean(dataCurr(indNoExt)); %Extra Northern
%             lpRtSoFlat = nanmean(dataCurr(indSoExt)); %Extra Southern
            
            if lpRtNoSteep > 0
                lpRtNoSteep = nan;
            end
            if lpRtSoSteep > 0
                lpRtSoSteep = nan;
            end
            if lpRtHiFlat > 0
                lpRtHiFlat = nan;
            end
            if lpRtNoFlat > 0
                lpRtNoFlat = nan;
            end
            if lpRtSoFlat > 0
                lpRtSoFlat = nan;
            end
            
            %Display lapse rates:
            disp([char(num2Month(zz)) ' lapse rates:']);
            disp(['     Northern, north facing, and steep = ' num2str(round2(1000*lpRtNoSteep,2)) ' C / km']);
            disp(['     Southern, south facing, and steep = ' num2str(round2(1000*lpRtSoSteep,2)) ' C / km']);
            disp(['     High and flat = '      num2str(round2(1000*lpRtHiFlat,2)) ' C / km']);
            disp(['     Northern and flat = ' num2str(round2(1000*lpRtNoFlat,2)) ' C / km']);
            disp(['     Southern and flat = ' num2str(round2(1000*lpRtSoFlat,2)) ' C / km']);
            
            lpRtCurr = nan(size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:))), 'single');

            lpRtCurr(indNoSteep) = lpRtNoSteep; %Northern and steep
            lpRtCurr(indSoSteep) = lpRtSoSteep; %Southern and steep
            lpRtCurr(indHiFlat ) = lpRtHiFlat; %High and flat
            lpRtCurr(indNoExt ) = lpRtNoFlat; %Extra Northern
            lpRtCurr(indSoExt ) = lpRtSoFlat; %Extra Southern

            
            %interpolate in order to fill in any nan's:
            if any2d(isnan(lpRtCurr))
                tempOut = lpRtCurr;
                
                indNan = find(isnan(tempOut));
                indData = setdiff(indAll, indNan);

                tempOut(indNan) = griddata(double(lonMeshLr(indData)), double(latMeshLr(indData)), ...
                    double(tempOut(indData)), ...
                    double(lonMeshLr(indNan)), double(latMeshLr(indNan)), sDs.intrp); 

                cntr = 0;
                while any2d(isnan(tempOut)) || cntr > 10
                    cntr = cntr + 1;
                    indNan = find(isnan(tempOut));
                    indData = setdiff(indAll, indNan);

                    tempOut(indNan) = griddata(double(lonMeshLr(indData)), double(latMeshLr(indData)), ...
                        double(tempOut(indData)), ...
                        double(lonMeshLr(indNan)), double(latMeshLr(indNan)), 'nearest');
                end
                
                lrLpsMnth(zz,:,:) = tempOut;
            else
                lrLpsMnth(zz,:,:) = lpRtCurr;
            end
        end
        clear zz


        %Smooth lapse rates using Gaussian filter (necessary to ensure no
        %discontinuities between regions):
        %NOTE: I do not fully understand the Gaussian filter code, but it seems to work!
        lrLpsMnthFilter = nan(size(lrLpsMnth), 'single');
        rad = 5;        
        sigma = rad/3; 

        xRad = (-rad : rad); 
        yRad = (-rad : rad)';  %+/-  1 grid cell (half each side)
        [xRadMesh, yRadMesh] = meshgrid(xRad,yRad);
        f = exp(-xRadMesh.^2/(2*sigma^2) - yRadMesh.^2/(2*sigma^2));
        f = f./sum(f(:));

        for zz = 1 : 12
            dataMean = mean2d(squeeze(lrLpsMnth(zz,:,:)));

            %Low pass filter:
            df = conv2(squeeze(lrLpsMnth(zz,:,:)) - dataMean, f, 'same');

            %Filtered data:
            lrLpsMnthFilter(zz,:,:) = df./conv2(ones(size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:)))),f,'same') + dataMean;
        end
        clear zz

        %Replace initial monthly lapse rates:
        lrLpsMnth = lrLpsMnthFilter;

        %Make figure of lapse rate zones:
        if blWrtFig == 1
            lpRtZones = nan(size(lpRtCurr));
            lpRtZones(indNoSteep) = 1; %Northern and steep
            lpRtZones(indSoSteep) = 2; %Southern and steep
            lpRtZones(indHiFlat ) = 3; %High and flat
            lpRtZones(indNoExt  ) = 4; %Extra Northern
            lpRtZones(indSoExt  ) = 5; %Extra Southern
            
            hFig = figure('color', 'white','visible','off'); hold on; 
            pcolor(sTInv{indLrDem}.(varLon), sTInv{indLrDem}.(varLat), double(lpRtZones)); colorbar; 
            hold off
            pathFig = fullfile(foldFig, [sDs.varDs '_lapse-rate-zones']);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
        %         export_fig([pathPlot '.eps'],'-painters');
            export_fig([pathFig '.png'],'-painters',strFigRes);
        end
    elseif strcmpi(lapseGrp, 'all')
        %%Calculate monthly lapse rates for all indices together:
        lrLpsMnth = nan(size(sTVar{sDs.indCm}.(varTLps)), 'single');
        for zz = 1 : 12
            dataCurr = squeeze(sTVar{sDs.indCm}.(varTLps)(zz,:,:));
            lpRtAll = nanmean(dataCurr(:));
            
            disp([char(num2Month(zz)) ' lapse rate = ' num2str(round2(1000*lpRtAll,2)) ' C / km']);
            
            lrLpsMnth(zz,:,:) = lpRtAll; %Northern and steep
        end
        clear zz
    else
        error('ERA_ds:unknownGrps',['The lapse rate group option ' strGrps ...
            ' has not been programmed for.']);
    end
    
    
    if blWrtFig == 1
        for mm = 1 : 12
            hFig = figure('color', 'white','visible','off'); hold on; 
            pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(squeeze(lrLpsMnth(mm,:,:)))); colorbar; 
            hold off
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [sDs.varDs '_lr_lapse-rate_' num2str(mm)]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
        %         export_fig([pathPlot '.eps'],'-painters');
            export_fig([pathFig '.png'],'-painters',strFigRes);
        end
        clear mm
    end
    
    
    %Interpolate smoothed low resolution lapse rate to high res:
    hrLpsMnth = nan([12, size(sTInv{indHrDem}.(varDem))], 'single');
    for mm = 1 : 12
        if regexpbl(sDs.intrp,'pchip')
                warning('off', 'PCHIP_2D:NaN');
            hrLpsMnth(mm,:,:) = ...
                PCHIP_2D(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), squeeze(lrLpsMnth(mm,:,:)), ...
                sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat));
                warning('on', 'PCHIP_2D:NaN');
        else
            hrLpsMnth(mm,:,:) = ...
                interp2(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), squeeze(lrLpsMnth(mm,:,:)), ...
                sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), sDs.intrp);
        end
    end
    clear mm



    %%Downscale (loop over unique month-year combinations, then process each day in that entry):
    %Define unique year-month combinations to loop over:
    datesYrMnth = unique(sTVar{sDs.indDs}.(varDate)(:,1:2), 'rows');
    %Keep only those specified for downscaling:
    datesYrMnth = datesYrMnth(datesYrMnth(:,1) >= min(sDs.yrsDs) & datesYrMnth(:,1) <= max(sDs.yrsDs), : );

    %Loop over output files:
    for ll = 1 : numel(datesYrMnth(:,1))
        %Current month:
        mnthWork = datesYrMnth(ll,2);
        yrWork = datesYrMnth(ll,1);

        nDy = eomday(datesYrMnth(ll,1), mnthWork);
        %Initialize grids for current month:
        indOutCurr = find(ismember(sTVar{sDs.indDs}.(varDate)(:,1:2), datesYrMnth(ll,1:2), 'rows'));
        
        sOutput = sTVar{sDs.indDs};
        sOutput.(sDs.varDs) = nan([nDy, size(sTInv{indHrDem}.(varDem))]);
        sOutput.('time') = sTVar{sDs.indDs}.('time')(indOutCurr);
        sOutput.(varDate) = sTVar{sDs.indDs}.(varDate)(indOutCurr,:);
        sOutput.(varLon) = sTInv{indHrDem}.(varLon);
        sOutput.(varLat) = sTInv{indHrDem}.(varLat);

        %Loop over days in month (for interpolation purposes):
        for dd = 1 : nDy
            %Find current date index:
            indDyCurr = find(ismember(sTVar{sDs.indDs}.(varDate), [datesYrMnth(ll,1:2), dd], 'rows'));

            %Only continue if data available for current date:
            if isempty(indDyCurr)
                continue
            end
            
            %Calculate daily weighting factor (used for lapse rate or
            %precipitation coefficient:
            indMnthPrv = mnthWork-1;
            if indMnthPrv == 0
                indMnthPrv = 12;
            end
            indMnthFut = mnthWork+1;
            if indMnthFut == 13
                indMnthFut = 1;
            end
            dyPrev = eomday(datesYrMnth(ll,1), indMnthPrv); %Editing year doesn't matter because Jan and Dec have same number of days every year.
            dyCurr = nDy/2;
            dyFut = eomday(datesYrMnth(ll,1), indMnthFut);
            xIn = [dyPrev/2, dyPrev + dyCurr/2, dyPrev + dyCurr + dyFut/2];
            xOut = dyPrev + dd;
            ind3Cons = [indMnthPrv, mnthWork, indMnthFut];

            %Interpolate to high-res DEM resolution: 
            dataDyCurrHr = interp2_norm_v2(...
                sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), sTInv{indLrDem}.(varArea), ...
                squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), ...
                sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varArea), ...
                sDs.intrp, strNorm);
            
            %Calculate a daily lapse rate grid (smoothed between three adjacent
            %months):
            %Interpolate between monthly lapse rates to get weighted daily values:
            hrLprDy = quad_fit_grids(squeeze(hrLpsMnth(ind3Cons,:,:)), xIn, xOut);

%                 warning('era_ds_himalaya:elevUnits','check elevation units for calculating lapse rate');
            %Apply lapse rate to interpolated temperature grid:
            sOutput.(sDs.varDs)(dd,:,:) = dataDyCurrHr + sTInv{indDemIntrp}.(varDemAnom).*hrLprDy;

            %Test and display conservation of mass:
            [sOutput.(sDs.varDs)(dd,:,:), wgtAvgIn, wgtAvgOut] =  norm_grid_v2(...
                squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), sTInv{indLrDem}.(varArea), ...
                sTVar{sDs.indDs}.(varLat), sTVar{sDs.indDs}.(varLon), ...
                squeeze(sOutput.(sDs.varDs)(dd,:,:)), sTInv{indHrDem}.(varArea), ...
                sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varLon), ...
                sDs.latDs, sDs.lonDs, strNorm);
%             [sOutput.(sDs.varDs)(dd,:,:), wgtAvgIn, wgtAvgOut] = norm_grid(squeeze(sTVar{sDs.indDs}.(varDs)(indDyCurr,:,:)), sTInv{indLrDem}.(varArea), ...
%                 sOutput.(sDs.varDs)(dd,:,:), sTInv{indHrDem}.(varArea), strNorm);

            if regexpbl(strNorm, 'mult')
                disp([num2str(yrWork) '-' num2str(mnthWork) '-' ...
                    num2str(dd) ': (Wgt Avg Out - Wgt Avg In) / Wgt Avg In = ' ...
                    num2str(round2(100*(wgtAvgOut-wgtAvgIn)/wgtAvgIn,1)) '%.']);
            elseif regexpbl(strNorm, 'add')
                disp([num2str(yrWork) '-' num2str(mnthWork) '-' ...
                    num2str(dd) ': Wgt Avg Out - Wgt Avg In = ' ...
                    num2str(round2(wgtAvgOut-wgtAvgIn,1)) '.']);
            elseif regexpbl(strNorm, 'none')
                disp([num2str(yrWork) '-' num2str(mnthWork) '-' ...
                    num2str(dd) ': bias = ' num2str(round2(wgtAvgOut-wgtAvgIn,1)) ' deg C']);
            else
                error('methodTlapse:unknownNorm',['The normalization method ' strNorm ' has not been prgorammed for.'])
            end


            if blWrtFig == 1 && dd == 15
                foldTran = fullfile(foldFig, 'transects');
                if ~exist(foldTran, 'dir')
                   mkdir(foldTran); 
                end
                
                %Make Transects
                rangeLon = max(sDs.lonDs) - min(sDs.lonDs);
                mnLonUse = min(sDs.lonDs) + 0.2*rangeLon;
                mxLonUse = max(sDs.lonDs) - 0.2*rangeLon;
                lonTrans = linspace(mnLonUse, mxLonUse, 3);
                for mm = 1 : numel(lonTrans)
                    [~, indLonHrDem] = min(abs(sTInv{indHrDem}.(varLon) - lonTrans(mm)));
                    [~, indLonLrDem] = min(abs(sTInv{indDemIntrp}.(varLon) - lonTrans(mm)));
                    [~, indLonSim] = min(abs(sTVar{sDs.indDs}.(varLon) - lonTrans(mm)));

                    transIn = squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,indLonSim));
                    transDsOut = squeeze(sOutput.(sDs.varDs)(dd,:,indLonHrDem));

                    hFig = figure('color', 'white', 'visible', 'off'); 
                    hold on

                    hTrans = plot(...
                        sTInv{   indHrDem}.(varLat), sTInv{   indHrDem}.(varDem)(:,indLonHrDem)/1000, '-',...
                        sTInv{indDemIntrp}.(varLat), sTInv{indDemIntrp}.(varDem)(:,indLonLrDem)/1000, '--',...
                        sTInv{indHrDem}.(varLat), transDsOut(:), '-', ...
                        sTInv{indLrDem}.(varLat), transIn(:), '-', ...
                        sTInv{indHrDem}.(varLat), dataDyCurrHr(:,indLonHrDem), '--', ...
                        'linewidth', 1.5);
                    cellLgd = {...
                        'Output Elevation (km)', ...
                        'Input Elevation, Interp (km)', ...
                        ['Output ' sDs.varDs ' (\circC)'], ...
                        ['Input ' sDs.varDs ' (\circC)'], ...
                        ['Input ' sDs.varDs ', Interp (\circC)'] ...
                        };
         
                    legend(hTrans, cellLgd, 'location', 'northwest');

                    ylabel(['Elevation (km) / ' sDs.varDs ' (\circC)']);
                    xlabel('Latitude (decimal degrees)')
                    hold off
                    pathFig = fullfile(foldTran, [sDs.varDs '_lon-' num2str(round2(lonTrans(mm),3)) ...
                        '_' num2str(yrWork) '-' num2str(mnthWork) '-' num2str(dd)]);
                    if exist([pathFig, '.fig'], 'file')
                        delete([pathFig, '.fig'])
                    end
                    if exist([pathFig, '.png'], 'file')
                        delete([pathFig, '.png'])
                    end
                    savefig(hFig, [pathFig '.fig']);
                    export_fig([pathFig '.png'],'-painters',strFigRes);
                end
                clear nn

            end
        end
        clear dd

    %     %Precip cannot be less than 0:
    %     if strcmpi(sDs.varDs, 'pr')
    %        output(output < 0) = 0; 
    %     end

        %Write monthly downscaled output to file:
        %Make month string:
        if mnthWork < 10
            strMnth = ['0' num2str(mnthWork)];
        else
            strMnth = num2str(mnthWork);
        end

        %Make file and path names:
        fileCurr = [sDs.varDs '_' sDs.timestep '_' strData '-ds-' sDs.region '_analysis_' ...
            num2str(datesYrMnth(ll,1)) strMnth '01' '-' ...
            num2str(datesYrMnth(ll,1)) strMnth num2str(nDy), '.nc'];
        
        ds_wrt_outputs(sOutput, 'output', sDs, sPath, 'file', fileCurr);
    end
    clear ll


    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes).' char(10)]);
end

warning('on', 'MATLAB:LargeImage');