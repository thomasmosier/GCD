function method_roe(sPath, sDs)

%Based on interpretation of methodology described in:
%Anders, A. M., Roe, G. H., Hallet, B., Montgomery, D. R., Finnegan, N. J.,
    %& Putkonen, J. (2006). Spatial patterns of precipitation and topography in
    %the Himalaya. Geological Society of America Special Papers, 398, 39?53.
%Roe, G. H. (2005). Orographic precipitation. Annu. Rev. Earth Planet.
    %Sci., 33, 645?671.



%Do not normalize grids:
strNorm = 'none';
% warning('methodTLapse:manualOpts','Set degN and elevHg in systematic/automatic way.')

wetThresh = 1.016; %Precipitation constituting a dry day (units = mm). Use 0 instead of USGS definition of 0.04 inches; https://earlywarning.usgs.gov/usraindry/rdreadme.php)

strRes = '300';


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
varTas = 'tas'; 
varVapSat = 'vappssat';

% varPW = 'tcw';
% varDeltaXPw = 'deltaxPW'; %Name to use for x change in PW
% varDeltaYPw = 'deltayPW'; %Name to use for y change in PW
% % varGradXPw = 'gradxPW'; %Name to use for x slope in PW
% % varGradYPw = 'gradyPW'; %Name to use for y slope in PW
% varDeltaMagPw = 'deltamagPW';

%General slopes/changes
varDeltaX = 'deltax'; %Name to use for x change
varDeltaY = 'deltay'; %Name to use for y change
varDeltaMag = 'deltamag'; %Name to use for combined gradient (sqrt[x^2+y^2])

%Gradients:
varUtmN = 'northing';
varUtmE = 'easting';
varUtm = 'utmzone';
varGradX = 'gradx'; %Name to use for x slope
varGradY = 'grady'; %Name to use for y slope
varGradMag = 'gradMag';
% varDeltaXAnom = 'deltaxAnom'; %Name to use for x slope
% varDeltaYAnom = 'deltayAnom'; %Name to use for y slope
% varDeltaMagAnom = 'deltaMagAnom';
varGradXAnom = 'gradxAnom'; %Name to use for x slope
varGradYAnom = 'gradyAnom'; %Name to use for y slope
varGradMagAnom = 'gradMagAnom';
% varGradYDemDiff = 'GradyElevAnom';
% varGradXDemDiff = 'GradxElevAnom';
% varDeltaMagDemDiff = 'GradMagElevAnom';
% varDemGradAnom = 'elevGradAnomMag';
varSlopeAnom = 'slopeAnom';

%Elevation
varDem = 'z'; %Data variable used in geopotential height field (also used to define elevation in other arrays)
varDemAnom = 'elevAnomMag'; %Variable to use for elevation anomaly
% varDeltaYDemAnom = 'deltayElevAnom';
% varDeltaXDemAnom = 'deltaxElevAnom';
varDeltaMagDemAnom = 'deltamagElevAnom';


if ~strcmpi(sDs.varDs, 'pr') && ~regexpbl(sDs.varDs, 'pre')
   error('methodTLapse:wrongVar',['The tlapse method is not designed to work with ' sDs.varDs]);
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


%Interpolate low resolution DEM and calculate gradients
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

%     %Calculate x,y changes in high-res elevation DEM anomaly (used for calculating alignment angle):
%     [sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom)] = gradient_map(sTInv{indDemIntrp}.(varDemAnom));
%     [~, sTInv{indDemIntrp}.(varDeltaMagDemAnom)] = cart2pol(sTInv{indDemIntrp}.(varDeltaXDemAnom), sTInv{indDemIntrp}.(varDeltaYDemAnom)); 
    

    [sTInv{indLrDem}.(varUtmN), sTInv{indLrDem}.(varUtmE), sTInv{indLrDem}.(varUtm)] = deg2utm_grid(sTInv{indLrDem}.(varLat), sTInv{indLrDem}.(varLon));

    [sTInv{indLrDem}.(varGradX), sTInv{indLrDem}.(varGradY), sTInv{indLrDem}.(varDeltaX), sTInv{indLrDem}.(varDeltaY), sTInv{indLrDem}.(varDeltaMag)] ...
        = gradient_utm_v3(sTInv{indLrDem}.(varDem), sTInv{indLrDem}.(varUtmN), sTInv{indLrDem}.(varUtmE), sTInv{indLrDem}.(varUtm));

    [sTInv{indHrDem}.(varUtmN), sTInv{indHrDem}.(varUtmE), sTInv{indHrDem}.(varUtm)] = deg2utm_grid(sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varLon));

    [sTInv{indHrDem}.(varGradX), sTInv{indHrDem}.(varGradY), sTInv{indHrDem}.(varDeltaX), sTInv{indHrDem}.(varDeltaY), sTInv{indHrDem}.(varDeltaMag)] ...
        = gradient_utm_v3(sTInv{indHrDem}.(varDem), sTInv{indHrDem}.(varUtmN), sTInv{indHrDem}.(varUtmE), sTInv{indHrDem}.(varUtm));

    szHrDem = size(sTInv{indHrDem}.(varDem));

    warning('off', 'gradientUtm:Nan');
    [sTInv{indDemIntrp}.(varGradX), sTInv{indDemIntrp}.(varGradY)] ...
        = gradient_utm_v3(sTInv{indDemIntrp}.(varDem), sTInv{indHrDem}.(varUtmN), sTInv{indHrDem}.(varUtmE), sTInv{indHrDem}.(varUtm));
    warning('on', 'gradientUtm:Nan');
    
    
    %Calculate elevation gradient magnitudes:
    [~,    sTInv{indLrDem}.(varGradMag)] = cart2pol(   sTInv{indLrDem}.(varGradX),    sTInv{indLrDem}.(varGradY));
%     [~, sTInv{indDemIntrp}.(varGradMag)] = cart2pol(sTInv{indDemIntrp}.(varGradX), sTInv{indDemIntrp}.(varGradY));
    [~,    sTInv{indHrDem}.(varGradMag)] = cart2pol(   sTInv{indHrDem}.(varGradX),    sTInv{indHrDem}.(varGradY));

%     sTInv{indDemIntrp}.(varDemGradAnom) = sTInv{indHrDem}.(varGradMag) - sTInv{indDemIntrp}.(varGradMag);

    %Calculate difference between interpolated DEM and HR DEM:

    warning('off', 'gradientUtm:Nan');
    [sTInv{indDemIntrp}.(varGradXAnom), sTInv{indDemIntrp}.(varGradYAnom)] ...
        = gradient_utm_v3(sTInv{indDemIntrp}.(varDemAnom), sTInv{indHrDem}.(varUtmN), sTInv{indHrDem}.(varUtmE), sTInv{indHrDem}.(varUtm));
    warning('on', 'gradientUtm:Nan');
    
%     [~, sTInv{indDemIntrp}.(varDeltaMagAnom)] = cart2pol(sTInv{indDemIntrp}.(varDeltaXAnom), sTInv{indDemIntrp}.(varDeltaYAnom));    
    [~, sTInv{indDemIntrp}.(varGradMagAnom)] = cart2pol(sTInv{indDemIntrp}.(varGradXAnom), sTInv{indDemIntrp}.(varGradYAnom));    

    sTInv{indDemIntrp}.(varSlopeAnom) = sTInv{indDemIntrp}.(varGradMagAnom) .* sign(sTInv{indDemIntrp}.(varDemAnom));
    
    
%     warning('off', 'gradientUtm:Nan');
%     [sTInv{indDemIntrp}.(varGradX), sTInv{indDemIntrp}.(varGradY), ~, ~, ~] ...
%         = gradient_utm_v3(sTInv{indDemIntrp}.(varDemAnom), sTInv{indHrDem}.(varUtmN), sTInv{indHrDem}.(varUtmE), sTInv{indHrDem}.(varUtm));
%     warning('on', 'gradientUtm:Nan');
%     [~, deltaMagDiffTemp] = cart2pol(sTInv{indDemIntrp}.(varGradX), sTInv{indDemIntrp}.(varGradY));
% 
%     sTInv{indDemIntrp}.(varDeltaMagDemDiff) = deltaMagDiffTemp;
%     indNeg = find(sTInv{indDemIntrp}.(varDemAnom) < 0);
%     sTInv{indDemIntrp}.(varDeltaMagDemDiff)(indNeg) = -sTInv{indDemIntrp}.(varDeltaMagDemDiff)(indNeg);
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
    export_fig([pathFig '.png'],'-painters',['-r' strRes]);
    
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
%     export_fig([pathFig '.png'],'-painters',['-r' strRes]);
    
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
    export_fig([pathFig '.png'],'-painters',['-r' strRes]);
    
%     %Magnitude of difference in SRTM and ERA DEMs
%     hFig = figure('color', 'white','visible','off'); hold on; 
%     histogram(sTInv{indDemIntrp}.(varDeltaMagDemAnom)(:), 100);
%     xlabel('Mag. Changes in High and Low DEMs (m)'); ylabel('Occurences');
%     pathFig = fullfile(foldFig, 'High-Low_res_DEM_anom_mag');
%     if exist([pathFig, '.fig'], 'file')
%         delete([pathFig, '.fig'])
%     end
%     if exist([pathFig, '.png'], 'file')
%         delete([pathFig, '.png'])
%     end
%     savefig(hFig, [pathFig '.fig']);
% %         export_fig([pathPlot '.eps'],'-painters');
%     export_fig([pathFig '.png'],'-painters',['-r' strRes]);
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
    
    
    %Calculate UTM coordinates:
    [sTVar{sDs.indDs}.(varUtmN), sTVar{sDs.indDs}.(varUtmE), sTVar{sDs.indDs}.(varUtm)] ...
        = deg2utm_grid(sTVar{sDs.indDs}.(varLat), sTVar{sDs.indDs}.(varLon));
    
    
    %CALCULATE SAPURATION VAPOR PRESSURE:
    indLrTas = nan;
    indHrTas = nan;
    for kk = 1 : numel(sDs.fldsTVar)
        if regexpbl(sDs.fldsTVar{kk}, {'tas','tave','tmp'}) && regexpbl(sDs.fldsTVar{kk}, {'sim','lr'})
            indLrTas = kk;
        elseif regexpbl(sDs.fldsTVar{kk}, {'tas','tave','tmp'}) && regexpbl(sDs.fldsTVar{kk}, {'out','hr'})
            indHrTas = kk;
        end
    end

    if isfield(sTVar{indLrTas}, ['att' varTas])
        unitsLrTas = find_att(sTVar{indLrTas}.(['att' varTas]), 'units');
        varTasLr = varTas;
    elseif isfield(sTVar{indLrTas}, 'atttmp')
        unitsLrTas = find_att(sTVar{indLrTas}.('atttmp'), 'units');
        varTasLr = 'tmp';
    else
        unitsLrTas = 'Celsius';
        varTasLr = varTas;
        warning('methodRoe:lrTasUnits',['The low resolution input ' ...
            'temperature units could not be found. They are being assumed to be degrees Celsius.']);
    end
    
    if isfield(sTVar{indHrTas}, ['att' varTas])
        unitsHrTas = find_att(sTVar{indHrTas}.(['att' varTas]), 'units');
        varTasHr = varTas;
    elseif isfield(sTVar{indHrTas}, 'atttmp')
        unitsHrTas = find_att(sTVar{indHrTas}.('atttmp'), 'units');
        varTasHr = 'tmp';
    else
        unitsHrTas = 'Celsius';
        varTasHr = varTas;
        warning('methodRoe:hrTasUnits',['The low resolution input ' ...
            'temperature units could not be found. They are being assumed to be degrees Celsius.']);
    end

    sTVar{indLrTas}.(varVapSat) = vapor_pressure_sat(sTVar{indLrTas}.(varTasLr), unitsLrTas);
    sTVar{indHrTas}.(varVapSat) = vapor_pressure_sat(sTVar{indHrTas}.(varTasHr), unitsHrTas);

    
    %Write inputs and bias corrected simulation:
    if ~isempty(sDs.wrtOut) && ~isempty(sDs.indRef) && ~isnan(sDs.indRef) 
%         ds_wrt_inputs(sTVar, sDs, sPath)
        ds_wrt_outputs(sTVar{sDs.indDs}, 'bc', sDs, sPath)
    end
    
    
    %Initialize low-res climatology structures:
    if ~isfield(sDs, 'indCm')
        sDs.indCm = numel(sTVar(:)) + 1;
    end
    sTVar{sDs.indCm} = struct;
        sTVar{sDs.indCm}.(varLat) = sTVar{sDs.indDs}.(varLat);
        sTVar{sDs.indCm}.(varLon) = sTVar{sDs.indDs}.(varLon);
        sTVar{sDs.indCm}.(varDate)= [nan([numel(mnthUse),1]), mnthUse(:), nan([numel(mnthUse),1])];
        sTVar{sDs.indCm}.('time') = nan([numel(mnthUse),1]);
        
        sTVar{sDs.indCm}.(sDs.varDs) = nan([numel(mnthUse), size(squeeze(sTVar{sDs.indDs}.(sDs.varDs)(1,:,:)))], 'single');
            
    %Calculate climatologies:
    wetDayGrid = zeros([numel(mnthUse), numel(sTVar{sDs.indDs}.(varLat)), numel(sTVar{sDs.indDs}.(varLon))], 'single');
    for mm = 1 : numel(mnthUse)
        %Precipitation climatology:
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
            export_fig([pathFig '.png'],'-painters',['-r' strRes]);
        end
        
%         %Precipitable water climatology:
%         %(wet days only):
%         for yy = 1 : numel(sTVar{sDs.indDs}.(varLat))
%             for xx = 1 : numel(sTVar{sDs.indDs}.(varLon))
%                 indWet = find(sTVar{sDs.indDs}.(sDs.varDs)(indCurr,yy,xx) >= wetThresh);
%                 
%                 if numel(indWet) > 1
%                     sTVar{sDs.indCm}.([sDs.varDs 'wet'])(mm,yy,xx) = nanmean(sTVar{sDs.indDs}.(sDs.varDs)(indCurr(indWet), yy, xx));
%                     wetDayGrid(mm,yy,xx) = numel(indWet);
%                 else
%                     sTVar{sDs.indCm}.(varPW)(mm,yy,xx) = 0; %valDryFlag;
%                 end
%             end
%         end
%         clear xx yy
        
        if blWrtFig == 1
            hFig = figure('color', 'white','visible','off'); hold on; 
            pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(squeeze(sTVar{sDs.indCm}.(sDs.varDs)(mm,:,:)))); colorbar; 
            hold off
            xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)');
            pathFig = fullfile(foldFig, [sDs.varDs '_wet-day_clim_' num2str(min(sDs.yrsBase)) '-' num2str(max(sDs.yrsBase)) '_' num2str(mnthUse(mm))]);
            if exist([pathFig, '.fig'], 'file')
                delete([pathFig, '.fig'])
            end
            if exist([pathFig, '.png'], 'file')
                delete([pathFig, '.png'])
            end
            savefig(hFig, [pathFig '.fig']);
            export_fig([pathFig '.png'],'-painters',['-r' strRes]);
        end
    end
    clear mm
    
 
    %Estimate orographic precipitation based on Roe et al.'s method
    %Fit parameters in Roe orographic precipitation model:
    wgtLagLr = neighbor_wgt_matrix(size(sTInv{indLrDem}.(varDem)), 'perp');
%         The first round of Roe orographic precipitation optimziation has occured. The MAE is 1.5342.
% This occurred for the paramter set = 0.0001, 0.0001, 0.012743, 1.6238.

    indYrs = find(sDs.yrsBase(1) <= sTVar{sDs.indDs}.(varDate)(:,1) & sTVar{sDs.indDs}.(varDate)(:,1) <= sDs.yrsBase(2));

    %Produce annual climatologies:
    annAvgSimVapSat = squeeze(nanmean(sTVar{indLrTas}.(varVapSat)(indYrs,:,:), 1)) / 100; %Vapor saturation in milibars
    annAvgSimPre = squeeze(nanmean(sTVar{sDs.indDs}.(sDs.varDs)(indYrs,:,:), 1));
%         annAvgDemGrad = squeeze(nanmean(sTVar{sDs.indDs}.(varWindDot), 1));

    nPts = 40; %Number of points used in parameter space for optimization

    prmOpt = roe_prm_opt(annAvgSimPre, annAvgSimVapSat, squeeze(sTInv{indLrDem}.(varGradMag)), wgtLagLr, nPts);

    [preModCov, preModOro] = roe_pre_empirical_v2(annAvgSimVapSat, squeeze(sTInv{indLrDem}.(varGradMag)), wgtLagLr, 'val', prmOpt);
        preModCov(preModCov < 0) = 0;
        preModOro(preModOro < 0) = 0;
    preModTot = preModCov + preModOro;


    if blWrtFig == 1
        %Roe convective precip:
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(preModCov)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Roe Convective Pre');
        pathFig = fullfile(foldFig, 'roe_convec_pre_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);

        %Roe orographic precip:
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(preModOro)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Roe Orographic Pre');
        pathFig = fullfile(foldFig, 'roe_orographic_pre_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);

        %Roe total pre:
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(preModTot)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Roe Total Pre');
        pathFig = fullfile(foldFig, 'roe_total_pre_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);
        
                %Fraction of orographic precip:
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(preModOro./preModTot)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Fraction of Roe Orographic Precipiation')
        pathFig = fullfile(foldFig, 'pre_model_frac_orographic_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);


        %Simulation input pre:
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(annAvgSimPre)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Simulation Total Pre');
        pathFig = fullfile(foldFig, 'sim_total_pre_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);

        %Modeled pre bias (Roe - input simulation pre):
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(preModTot-annAvgSimPre)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Precip Model Bias (Roe - Sim)');
        pathFig = fullfile(foldFig, 'pre_model_bias_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);

        %Percent Error in Modeled pre (Roe - input simulation pre):
        hFig = figure('color', 'white','visible','off'); hold on; 
        pcolor(sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), double(100*(preModTot-annAvgSimPre)./annAvgSimPre)); colorbar; 
        hold off;
        xlabel('Longitude (degrees East)'); ylabel('Latitude (degrees North)'); title('Precip Model Percent Error ([Roe - Sim] / Sim)');
        pathFig = fullfile(foldFig, 'pre_model_percent_error_ann_clim');
        if exist([pathFig, '.fig'], 'file')
            delete([pathFig, '.fig']);
        end
        if exist([pathFig, '.png'], 'file')
            delete([pathFig, '.png']);
        end
        savefig(hFig, [pathFig '.fig']);
    %         export_fig([pathPlot '.eps'],'-painters');
        export_fig([pathFig '.png'],'-painters',['-r' strRes]);
    end  

    
    %Calculate high-res Weight matrix for downscaling:
    wgtLagHr = neighbor_wgt_matrix(size(sTInv{indHrDem}.(varDem)), 'perp');
    
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
           
            %Interpolate precip to output resolution:
            dataDyCurrHr = interp2_norm_v2(...
                sTVar{sDs.indDs}.(varLon), sTVar{sDs.indDs}.(varLat), sTInv{indLrDem}.(varArea), ...
                squeeze(sTVar{sDs.indDs}.(sDs.varDs)(indDyCurr,:,:)), ...
                sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), sTInv{indHrDem}.(varArea), ...
                sDs.intrp, strNorm);
            dataDyCurrHr(dataDyCurrHr < 0) = 0;

            %Get saturation vapor pressure for current day (previously based on downscaled
            %temperature)
            if min(sTInv{indHrDem}.(varLat)) > min(sTVar{indHrTas}.(varLat)) || max(sTInv{indHrDem}.(varLat)) < max(sTVar{indHrTas}.(varLat)) 
                error('methodRoe:hrTasLatOutBnds', ['The high-resolution '...
                    'temperature latitude grid is larger than the input '...
                    'high-resolution DEM. This has not been programmed for.']);
            elseif min(sTInv{indHrDem}.(varLon)) > min(sTVar{indHrTas}.(varLon)) || max(sTInv{indHrDem}.(varLon)) < max(sTVar{indHrTas}.(varLon)) 
                error('methodRoe:hrTasLonOutBnds', ['The high-resolution '...
                    'temperature longitude grid is larger than the input '...
                    'high-resolution DEM. This has not been programmed for.']);
            end
            
            vaporSatCurr = nan(szHrDem,'single');
            
            rowVapUse = find(ismember(round2(sTInv{indHrDem}.(varLat), 4), round2(sTVar{indHrTas}.(varLat), 4)) == 1);
            colVapUse = find(ismember(round2(sTInv{indHrDem}.(varLon), 4), round2(sTVar{indHrTas}.(varLon), 4)) == 1);
                [colVapUseMesh, rowVapUseMesh] = meshgrid(colVapUse(:)', rowVapUse(:));
                indVapUse = sub2ind(szHrDem, rowVapUseMesh(:), colVapUseMesh(:));
            vaporSatCurr(indVapUse) = squeeze(sTVar{indHrTas}.(varVapSat)(dd,:,:));

            
            %Estimate orographic precipitation based on unaccounted for
            %slope:
            [~, oroAnomHr] = roe_pre_empirical_v2(vaporSatCurr, sTInv{indDemIntrp}.(varSlopeAnom), wgtLagHr, 'val', prmOpt);

            %Estimate orographic precip at high-resolution:
            [conDyHr, oroDyHr] = roe_pre_empirical_v2(vaporSatCurr, sTInv{indHrDem}.(varGradMag), wgtLagHr, 'val', prmOpt);
            
            %Unaccounted for precipitation fraction:
            oroAnomFact = oroAnomHr ./ (oroDyHr + conDyHr);
                oroAnomFact(oroDyHr + conDyHr == 0) = 0;

            %Add orographic precip factor to interpolated precip:
            dataOutTemp = dataDyCurrHr.*(1 + oroAnomFact);
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
            

            %Write test figure for each month of first year
            if blWrtFig == 1 && dd == 15 && yrWork == datesYrMnth(1,1)
                %Create vector of output times:
                timeVecNan = nan(1,3);

                fileOroPr = ['orographic_pr_anom_factor_' reg '_' num2str(yrWork) '-' num2str(mnthCurr) '-' num2str(dd) '.nc'];
                pathOroPr = fullfile(foldFig, fileOroPr);
                %If file already exists, delete:
                if exist(pathOroPr, 'file')
                    delete(pathOroPr);
                end

                warning('off', 'write_NC:nanTime');
                print_grid_NC_v2(pathOroPr, oroAnomFact, 'orographic_pr_anom_factor', sTInv{indHrDem}.(varLon), sTInv{indHrDem}.(varLat), ...
                    timeVecNan, timeVecNan);
                warning('on', 'write_NC:nanTime');
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
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes).' char(10)]);
end

warning('on', 'MATLAB:LargeImage');