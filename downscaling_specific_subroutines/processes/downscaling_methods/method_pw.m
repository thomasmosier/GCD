function method_pw(sPath, sDs)

%Do not normalize grids:
strNorm = 'none';
% warning('methodTLapse:manualOpts','Set degN and elevHg in systematic/automatic way.')

wetThresh = 1.016; %Precipitation constituting a dry day (units = mm). Use 0 instead of USGS definition of 0.04 inches; https://earlywarning.usgs.gov/usraindry/rdreadme.php)

strFigRes = '-r600';

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

%Precipitable water variables
varPW = 'tcw';
varDeltaXPw = 'deltaxPW'; %Name to use for x change in PW
varDeltaYPw = 'deltayPW'; %Name to use for y change in PW
% varGradXPw = 'gradxPW'; %Name to use for x slope in PW
% varGradYPw = 'gradyPW'; %Name to use for y slope in PW
varDeltaMagPw = 'deltamagPW';

%General slopes/changes
varDeltaX = 'deltax'; %Name to use for x change
varDeltaY = 'deltay'; %Name to use for y change
varDeltaMag = 'deltamag'; %Name to use for combined gradient (sqrt[x^2+y^2])

%Elevation
varDem = 'z'; %Data variable used in geopotential height field (also used to define elevation in other arrays)
varDemAnom = 'elevAnomMag'; %Variable to use for elevation anomaly
varDeltaYDemAnom = 'deltayElevAnom';
varDeltaXDemAnom = 'deltaxElevAnom';
varDeltaMagDemAnom = 'deltamagElevAnom';

if ~strcmpi(sDs.varDs, 'pr') && ~regexpbl(sDs.varDs, 'pre')
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
sTInv = ds_ld_fields(sPath, sDs.fldsTInv, sDs.lonDs, sDs.latDs, nan(1,2), nan(1,1), 'onefile');

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
        sTInv{indDemIntrp}.(varDem) = interp2(...
            lonMeshLr, latMeshLr, ...
            squeeze(sTInv{indLrDem}.(varDem)), ...
            lonMeshHr, latMeshHr, sDs.intrp);
    end

    %Calculate difference between interpolated DEM and HR DEM:
    %This is used later on for determining lapse rate correct to apply
    sTInv{indDemIntrp}.(varDemAnom) = sTInv{indHrDem}.(varDem) - sTInv{indDemIntrp}.(varDem);

    %Calculate x,y changes in high-res elevation DEM anomaly (used for calculating alignment angle):
    [sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom)] = gradient_map(sTInv{indDemIntrp}.(varDemAnom));
    [~, sTInv{indDemIntrp}.(varDeltaMagDemAnom)] = cart2pol(sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom));    
else
    error('methodTlapse:DemInputsNotFound','Either the low or high resolution DEM inputs have not been located.');
end


%Plot elevation histograms
if blWrtFig == 1
    %Low-res elevation change magnitude
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
       error('methodPw:missingMonths', ['The method_pw function is '...
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
    
    %Find which time-series input is precipitable water (pw):
    if numel(sTVar(:)) > 0
        indPw = [];
        for jj = 1 : numel(sDs.fldsTVar)
           if regexpbl(sDs.fldsTVar{jj}, {'pw', varPW}) 
               indPw(end+1) = jj;
           end
        end
        clear jj
        
        if numel(indPw) ~= 1
           error('methodPW:noPwVar',[numel(indPw) ' was found, but there should be exactly 1.']) 
        end
        
        %Check that reference coordinates are equal and add PW to main sTVar{sDs.indDs} array:
        if ~isequal(size(sTVar{indPw}.(varPW)), size(sTVar{sDs.indDs}.(sDs.varDs))) || ~isequal(sTVar{indPw}.(varLon), sTVar{sDs.indDs}.(varLon)) || ~isequal(sTVar{indPw}.(varLat), sTVar{sDs.indDs}.(varLat)) || ~isequal(sTVar{indPw}.time, sTVar{sDs.indDs}.time)
            error('methodPw:TVarDiffSize','The input precipitable water and precipitation time-series have different sizes.');
        end
    end
    
    
    %Write inputs and bias corrected simulation:
    if ~isempty(sDs.wrtOut) && ~isempty(sDs.indRef) && ~isnan(sDs.indRef) 
%         ds_wrt_inputs(sTVar, sDs, sPath)
        ds_wrt_outputs(sTVar{sDs.indDs}, 'bc', sDs, sPath)
    end
    
    
    %Initialize climatology structures:
    %Low-res:
    if ~isfield(sDs, 'indCm')
        sDs.indCm = numel(sTVar(:)) + 1;
    end
    sTVar{sDs.indCm} = struct;
        sTVar{sDs.indCm}.(varLat) = sTVar{sDs.indDs}.(varLat);
        sTVar{sDs.indCm}.(varLon) = sTVar{sDs.indDs}.(varLon);
        sTVar{sDs.indCm}.(varDate)= [nan([numel(mnthUse),1]), mnthUse(:), nan([numel(mnthUse),1])];
        sTVar{sDs.indCm}.('time') = nan([numel(mnthUse),1]);
        
        sTVar{sDs.indCm}.(sDs.varDs) = nan([numel(mnthUse), size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:)))], 'single');
        sTVar{sDs.indCm}.(varDeltaXPw)   = sTVar{sDs.indCm}.(sDs.varDs);
        sTVar{sDs.indCm}.(varDeltaYPw)   = sTVar{sDs.indCm}.(sDs.varDs);
        sTVar{sDs.indCm}.(varDeltaMagPw) = sTVar{sDs.indCm}.(sDs.varDs);
        
    %Interpolated:
    if ~isfield(sDs, 'indCmIntrp')
        sDs.indCmIntrp = numel(sTVar(:)) + 1;
    end
    sTVar{sDs.indCmIntrp} = struct;
        sTVar{sDs.indCmIntrp}.(varLat) = sTInv{indHrDem}.(varLat);
        sTVar{sDs.indCmIntrp}.(varLon) = sTInv{indHrDem}.(varLon);
        sTVar{sDs.indCmIntrp}.(varDate) = [nan([numel(mnthUse),1]), mnthUse(:), nan([numel(mnthUse),1])];
        sTVar{sDs.indCmIntrp}.('time') = nan([numel(mnthUse),1]);
        
        sTVar{sDs.indCmIntrp}.(varDeltaXPw) = nan([numel(mnthUse), size(sTInv{indHrDem}.(varDem))], 'single');
        sTVar{sDs.indCmIntrp}.(varDeltaYPw)        = sTVar{sDs.indCmIntrp}.(varDeltaXPw);
        sTVar{sDs.indCmIntrp}.(varDeltaMagPw)      = sTVar{sDs.indCmIntrp}.(varDeltaXPw);
        sTVar{sDs.indCmIntrp}.(varDeltaMagDemAnom) = sTVar{sDs.indCmIntrp}.(varDeltaXPw);
        sTVar{sDs.indCmIntrp}.(varDemAnom)         = sTVar{sDs.indCmIntrp}.(varDeltaXPw);
            
    %Calculate climatologies:
    wetDayGrid = zeros([numel(mnthUse), numel(sTVar{indPw}.(varLat)), numel(sTVar{indPw}.(varLon))], 'single');
    for mm = 1 : numel(mnthUse)
        %Precipitation climatology:
        indCurr = find(sTVar{sDs.indDs}.(varDate)(:,1) >= min(sDs.yrsBase) & sTVar{sDs.indDs}.(varDate)(:,1) <= max(sDs.yrsBase) & sTVar{sDs.indDs}.(varDate)(:,2) == mnthUse(mm));
        sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:) = squeeze(nanmean(sTVar{sDs.indDs}.(sDs.varDs)(indCurr,:,:), 1));
        
        if blWrtFig == 1
            hFig = figure('color', 'white','visible','off'); hold on; 
            hPClr = imagesc(sTVar{sDs.indCm}.(varLon), sTVar{sDs.indCm}.(varLat), double(squeeze(sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:)))); colorbar; 
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
        
        %Precipitable water climatology:
        %(wet days only):
        for yy = 1 : numel(sTVar{indPw}.(varLat))
            for xx = 1 : numel(sTVar{indPw}.(varLon))
                indWet = find(sTVar{sDs.indDs}.(sDs.varDs)(indCurr,yy,xx) >= wetThresh);
                
                if numel(indWet) > 1
                    sTVar{sDs.indCm}.(varPW)(mm,yy,xx) = nanmean(sTVar{indPw}.(varPW)(indCurr(indWet), yy, xx));
                    wetDayGrid(mm,yy,xx) = numel(indWet);
                else
                    sTVar{sDs.indCm}.(varPW)(mm,yy,xx) = 0; %valDryFlag;
                end
            end
        end
        clear xx yy
        
        if blWrtFig == 1
            hFig = figure('color', 'white','visible','off'); hold on; 
            imagesc(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(squeeze(sTVar{sDs.indCm}.(varPW)(mm,:,:)))); colorbar; 
            hold off
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [varPW '_wet-day_clim_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
            export_fig([pathFig '.png'],'-painters',strFigRes);
        end
    end
    clear mm
    
    
    %Calculate other climatological parameters (e.g. spatial changes):
    for mm = 1 : numel(mnthUse)
        [sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:), sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)] ...
            = gradient_map(squeeze(sTVar{sDs.indCm}.(varPW)(mm,:,:)));
        
%         %Set gradients for grid cells where there are no wet days to nan:
%         sTVar{sDs.indCm}.(varDeltaXPw)(isnan(sTVar{sDs.indCm}.(varPW))) = nan;
%         sTVar{sDs.indCm}.(varDeltaYPw)(isnan(sTVar{sDs.indCm}.(varPW))) = nan;

        qLr = -alignment_factor(...
            sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY), ...
            squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)), squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)));

        tempDeltaXLr = squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)) ./ squeeze(sTVar{sDs.indCm}.(varPW)(mm,:,:));
        tempDeltaYLr = squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)) ./ squeeze(sTVar{sDs.indCm}.(varPW)(mm,:,:));
            tempDeltaXLr(isnan(tempDeltaXLr) & ~isnan(squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)))) = 0;
            tempDeltaYLr(isnan(tempDeltaYLr) & ~isnan(squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)))) = 0;
        sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:) = tempDeltaXLr;
        sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:) = tempDeltaYLr;
        
        %Scale combined climatological change by alignment factor, q.
        [~, tempMagOut] = cart2pol(...
            squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)), squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)));
        sTVar{sDs.indCm}.(varDeltaMagPw)(mm,:,:) = qLr.*tempMagOut;
    
        
        %If precipitation, calculate high-res gradient as well:
        %Interpolate PW changes and gradients to high-res grid:
        tempDeltaXHr = interp2(...
            sTVar{sDs.indCm}.(varLon), sTVar{sDs.indCm}.(varLat), ...
            squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)), ...
            sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), sDs.intrp);
        tempDeltaYHr = interp2(...
            sTVar{sDs.indCm}.(varLon), sTVar{sDs.indCm}.(varLat), ...
            squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)), sTInv{indHrDem}.(varLon), ...
            sTInv{indHrDem}.(varLat), sDs.intrp);
%         tempDeltaXHr(isnan(tempDeltaXHr)) = 0;
%         tempDeltaYHr(isnan(tempDeltaYHr)) = 0;
        
        sTVar{sDs.indCmIntrp}.(varDeltaXPw)(mm,:,:) = tempDeltaXHr;
        sTVar{sDs.indCmIntrp}.(varDeltaYPw)(mm,:,:) = tempDeltaYHr;
        

%         %Calculate magnitude of change for interpolated PW
%         [~, sTVar{sDs.indCmIntrp}.(varDeltaMagPw)(mm,:,:)] = cart2pol(...
%             squeeze(sTVar{sDs.indCmIntrp}.(varDeltaXPw)(mm,:,:)), ...
%             squeeze(sTVar{sDs.indCmIntrp}.(varDeltaYPw)(mm,:,:)));

        %Calculate alignment between elevation anomaly gradients and
        %interpolated PW gradients
        %For precip, switch sign because decreasing PW corresponds to
        %increasing precip
        qHr = -alignment_factor(sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom), ...
            squeeze(sTVar{sDs.indCmIntrp}.(varDeltaXPw)(mm,:,:)), squeeze(sTVar{sDs.indCmIntrp}.(varDeltaYPw)(mm,:,:)));

        %Calculate elevation difference in direction of interpolated PW
        %fractional gradient:
%         sTVar{sDs.indCmIntrp}.(varDemAnom)(mm,:,:) = qHr.*sTInv{indDemIntrp}.(varDemAnom);
        sTVar{sDs.indCmIntrp}.(varDeltaMagDemAnom)(mm,:,:) = qHr.*sTInv{indDemIntrp}.(varDeltaMagDemAnom);
        
        %Produce figure to check processing:
        if blWrtFig == 1
            [lrPrGridLon,   lrPrGridLat] = meshgrid(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat));
            [lrDemGridLon, lrDemGridLat] = meshgrid(sTInv{ indLrDem}.(varLon), sTInv{ indLrDem}.(varLat));

            %PW gradient with arrows:
            hFig = figure('color', 'white','visible','off'); 
            hold on; 
            imagesc(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(qLr)); colorbar; 
            quiver(lrPrGridLon, lrPrGridLat, ...
                squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)), ...
                squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)),...
                'color', 'black'); 
            quiver(lrDemGridLon, lrDemGridLat, ...
                sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY), 'color', 'red');
            text(double(0.5*mean(diff(lrPrGridLon(1,:)))+lrPrGridLon(2,2)), ...
                double(0.5*mean(diff(lrPrGridLon(:,1)))+lrPrGridLat(2,2)), ...
                'Black = PW changes; Red = Elevation changes', 'FontSize',14);
            hold off;
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [varPW '-DEM_grad_alignment_' num2str(mm)]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
        %         export_fig([pathPlot '.eps'],'-painters');
            export_fig([pathFig '.png'],'-painters',strFigRes);

            %Alignment with arrows:
            hFig = figure('color', 'white','visible','off'); hold on; 
            imagesc(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), ...
                double(squeeze(sTVar{sDs.indCm}.(varDeltaMagPw)(mm,:,:)))); colorbar; 
            quiver(lrPrGridLon, lrPrGridLat, ...
                squeeze(sTVar{sDs.indCm}.(varDeltaXPw)(mm,:,:)), ...
                squeeze(sTVar{sDs.indCm}.(varDeltaYPw)(mm,:,:)), ...
                'color', 'black'); 
            quiver(lrDemGridLon, lrDemGridLat, ...
                sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY), 'color', 'red');
            text(double(0.5*mean(diff(lrPrGridLon(1,:)))+lrPrGridLon(2,2)), ...
                double(0.5*mean(diff(lrPrGridLon(:,1)))+lrPrGridLat(2,2)), ...
                'Black = PW changes; Red = Elevation changes', 'FontSize',14);
            hold off;
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [varPW '_grad_magnitude_' num2str(mm)]);
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
    
    
    %Synthesize regression equations that are the basis of downscaling:
    cfReg = nan(12,2);
    rSq   = nan(12,1);
    disp('Calculating empirical relationship between changes in precipitable water and elevation change.');
    disp('Units = fractional change in precipitable water per km change in elevation');
    for mm = 1 : numel(mnthUse)
        %Current climate variable to be predicted = magnitude of the
        %fractional PW gradients in the direction of the ERA GP gradient
        pwDeltaCurr = squeeze(sTVar{sDs.indCm}.(varDeltaMagPw)(mm,:,:));
        
        %Only use data points where there are wet days
        pwDeltaCurr(squeeze(wetDayGrid(mm,:,:)) == 0) = nan;

        indReg = find(~isnan(sTInv{indLrDem}.(varDeltaMag)(:)) & ~isnan(pwDeltaCurr(:)));

        %Fit linear function (input order is x, y):
        cfReg(mm,:) = polyfit(sTInv{indLrDem}.(varDeltaMag)(indReg), pwDeltaCurr(indReg), 1);
        rSq(mm)     = pearson_cf(cfReg(mm, :), sTInv{indLrDem}.(varDeltaMag)(indReg), pwDeltaCurr(indReg));

        %Can alternatively use multiple linear regression (edited version of
        %function from File Exchange). Answers somewhat similar.
        %Do linear regression (input order is y, x; can have multiple columns for x):

%         [cfRegm(kk), RSq(kk)] ...
%           = mregress_condition(areaLr(indReg).*climCurr(indReg), areaLr(indReg).*sTInv{indLrDem}.(varGradMag)(indReg), 0);

        if round2(cfReg(mm, 2),1) ~= 0.0
            warning('ERA_ds:yIntNonZero', ...
                ['The intertercept has been estimated as ' ...
                num2str(cfReg(mm, 2)) '. It should be zero.']);
        end
        
        disp(['month ' num2str(mnthUse(mm)) ...
            ' (using ' num2str(numel(indReg)) ' data points): coef = ' ...
            num2str(round2(1000*cfReg(mm), 2)) '; r^2 = ' num2str(round2(rSq(mm), 2))]);


        %Make figures of corrections:
        if blWrtFig == 1
            %Scatter plot of regression inputs:
            hFig = figure('color', 'white','visible','off'); hold on; 
            scatter(sTInv{indLrDem}.(varDeltaMag)(indReg), pwDeltaCurr(indReg));
            xlabel('Change in Elevation (km)'); ylabel('Change in Total Column Water (fraction)');
            pathFig = fullfile(foldFig, ['tcw_elevation_regression_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
            export_fig([pathFig '.png'],'-painters',strFigRes);


            %Low-res histogram
            hFig = figure('color', 'white','visible','off'); hold on; 
            histogram(pwDeltaCurr(indReg));
            xlabel('Magnitude of Changes in PW'); ylabel('Occurences');
            pathFig = fullfile(foldFig, [varPW '_fractional_grad_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
            export_fig([pathFig '.png'],'-painters',strFigRes);

            %High-res histogram
            tempHrCorr = squeeze(sTVar{sDs.indCmIntrp}.(varDeltaMagDemAnom)(mm,:,:));
            hFig = figure('color', 'white','visible','off'); hold on; 
            histogram(cfReg(mm,1)*tempHrCorr(:),100);
            xlabel('PW Empirical Correction'); ylabel('Occurences');
            pathFig = fullfile(foldFig, [varPW '_hr_corrections_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
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
            
            %Smoothed regression factor for current day:
            cfPreDy = quad_fit_grids(cfReg(ind3Cons,1), xIn, xOut);

            %Interpolate precip to output resolution:
            [dataDyCurrHr, ~, ~] = interp2_norm_v2( ...
                sTInv{indLrDem}.(varLon), sTInv{indLrDem}.(varLat), sTInv{indLrDem}.(varArea),  ...
                squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), ...
                sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varArea), ...
                sDs.intrp, strNorm);
            %Values cannot be negative
            dataDyCurrHr(dataDyCurrHr < 0) = 0;

            %Interpolate monthly PW gradient and coefficients to daily:
            hrZPWGradDy = quad_fit_grids(sTVar{sDs.indCmIntrp}.(varDeltaMagDemAnom)(ind3Cons,:,:), xIn, xOut); %This is the elevation anomaly gradient in the direction of the PW climatology gradient

            %Estimate daily precipitation values:
            dataOutTemp = dataDyCurrHr.*(1 + cfPreDy*squeeze(hrZPWGradDy));
            %Values cannot be negative
            dataOutTemp(dataOutTemp < 0)= 0;
             
            %Assign to output array:
            sOutput.(sDs.varDs)(dd,:,:) = dataOutTemp;

            %Test and display conservation of mass:
            [sOutput.(sDs.varDs)(dd,:,:), wgtAvgIn, wgtAvgOut] =  norm_grid_v2(...
                squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), sTInv{indLrDem}.(varArea), ...
                sTVar{sDs.indDs}.(varLat), sTVar{sDs.indDs}.(varLon), ...
                squeeze(sOutput.(sDs.varDs)(dd,:,:)), sTInv{indHrDem}.(varArea), ...
                sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varLon), ...
                sDs.latDs, sDs.lonDs, strNorm);
%             [sOutput.(sDs.varDs)(dd,:,:), wgtAvgIn, wgtAvgOut] = norm_grid(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), sTInv{indLrDem}.(varArea), ...
%                 sOutput.(sDs.varDs)(dd,:,:), sTInv{indHrDem}.(varArea), strNorm);
            if ~isempty(wgtAvgIn) && ~isempty(wgtAvgOut)
                if regexpbl(strNorm, 'mult')
                    disp([num2str(yrWork) '-' num2str(mnthWork) '-' ...
                        num2str(dd) ': (Wgt Avg Out - Wgt Avg In) / Wgt Avg In = ' ...
                        num2str(round2(100*(wgtAvgOut-wgtAvgIn)/wgtAvgIn,1)) '%.']);
                elseif regexpbl(strNorm, 'none')
                    disp([num2str(yrWork) '-' num2str(mnthWork) '-' ...
                        num2str(dd) ': ratio = ' num2str(round2(wgtAvgOut/wgtAvgIn,2)) ]);
                else
                    error('methodTlapse:unknownNorm',['The normalization method ' strNorm ' has not been prgorammed for.'])
                end
            else
                disp([num2str(yrWork) '-' num2str(mnthWork) ' processed.']);
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
                    
                    %Find scaling factor
                    tranScl = max([round2(max(sTInv{indDemIntrp}.(varDem)(:,indLonLrDem)/1000) / max(transDsOut(:)), -1),1]);
                    
                    %Make plot
                    hTrans = plot(...
                        sTInv{   indHrDem}.(varLat), sTInv{   indHrDem}.(varDem)(:,indLonHrDem)/1000, '-',...
                        sTInv{indDemIntrp}.(varLat), sTInv{indDemIntrp}.(varDem)(:,indLonLrDem)/1000, '--',...
                        sTInv{indHrDem}.(varLat), tranScl*transDsOut(:), '-', ...
                        sTInv{indLrDem}.(varLat), tranScl*   transIn(:), '-', ...
                        sTInv{indHrDem}.(varLat), tranScl*dataDyCurrHr(:,indLonHrDem), '--', ...
                        'linewidth', 1.5);
                    cellLgd = {...
                        'Output Elevation (km)', ...
                        'Input Elevation, Interp (km)', ...
                        ['Output ' sDs.varDs ' (' num2str(tranScl) '*mm)'], ...
                        ['Input ' sDs.varDs ' (' num2str(tranScl) '*mm)'], ...
                        ['Input, ' sDs.varDs ', Interp (' num2str(tranScl) '*mm)'], ...
                        };
%                     else
%                         hTrans(end+1) = plot(sTVar{sDs.indDs}.(varLat), tranScl*transIn(:), 'linewidth', 1.5);
%                         cellLgd{end+1} = ['ERA ' sDs.varDs ' (' num2str(tranScl) '*' unitsOut ')'];
                    legend(hTrans, cellLgd, 'location', 'northwest');

                    ylabel(['Elevation (km) / ' sDs.varDs ' (' num2str(tranScl) '*mm)']);
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
            num2str(datesYrMnth(ll,1)) strMnth num2str(nDy) '.nc'];
        
        ds_wrt_outputs(sOutput, 'output', sDs, sPath, 'file', fileCurr);
    end
    clear ll


    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / sDs.nLp;
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes).' char(10)]);
end

warning('on', 'MATLAB:LargeImage');