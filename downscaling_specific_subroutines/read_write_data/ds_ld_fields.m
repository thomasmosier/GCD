function [sData, varargout] = ds_ld_fields(sPath, fldsLd, lonDs, latDs, yrsLd, mnthsLd, varargin)

unitsTable = cell(0,2);
blStitch = 1;
indiceIn = nan;
resample = 'none';
fr = 4;
cropType = 'out';
filesLd = '';
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'units')
            unitsTable = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'stitch')
            blStitch = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'indice')
            indiceIn = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'resample')
            resample = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'frame')
            fr = varargin{ii+1};
        elseif strcmpi(varargin{ii}, 'onefile')
            filesLd = 'onefile';
        end
    end
end

varName = 'description';
varLat = 'latitude';
varLon = 'longitude';
varTime = 'time';
varDate = 'date';


for ii = numel(fldsLd) : -1 : 1
    if regexpbl(fldsLd{ii}, 'bc') && ~regexpbl(fldsLd{ii}, '4bc') 
        fldsLd(ii) = [];
    end
end


sData = cell(numel(fldsLd), 1);
indSimProj = [];
indSimHist = [];
for ii = 1 : numel(fldsLd)
    indUnd = regexpi(fldsLd{ii}, '_');
    if ~isempty(indUnd)
        varLd = fldsLd{ii}(indUnd(end)+1:end);
    else
        warning('methodDirect:unknownTimeInvFld', [char(39) fldsLd{ii} char(39) ' has not been programmed for.']);
    end

    if all(~isnan(mnthsLd))
        sData{ii} = read_geodata_v2(char(sPath.(fldsLd{ii})), varLd, ...
            lonDs, latDs, yrsLd, fr, cropType, 'months', mnthsLd, filesLd);
    else
        sData{ii} = read_geodata_v2(char(sPath.(fldsLd{ii})), varLd, ...
            lonDs, latDs, yrsLd, fr, cropType, filesLd);
    end
    sData{ii}.(varName) = fldsLd{ii};
    
    %Check that data for all years loaded:
    if isfield(sData{ii}, varDate)
        yrsPresent = unique(sData{ii}.(varDate)(:,1));
        if ~all(isnan(yrsPresent)) && ~isempty(setdiff((min(yrsLd):max(yrsLd)), yrsPresent)) && max(yrsLd) ~= max(yrsPresent) && min(yrsLd) ~= min(yrsPresent) && any(diff(yrsPresent) ~= 1)
           warning('dsLdFields:diffNumberYrs', ['The number of years loaded '...
               'does not equal the number requested for ' sPath.(fldsLd{ii}) '.']); 
        end
    end
    
    warning('off', 'MATLAB:mode:EmptyInput')
    if isfield(sData{ii}, 'date')
        if all(ismember(sData{ii}.('date')(:,2), mnthsLd) ~= 0) && (numel(sData{ii}.('date')(1,:)) == 2 || all(sData{ii}.('date')(:,3) == sData{ii}.('date')(1,3)))
            sData{ii}.timestep = 'monthly';
        elseif numel(sData{ii}.('date')(1,:)) == 3 && mode(abs(diff(sData{ii}.date(:,3)))) == 1
            sData{ii}.timestep = 'daily';
        else
            sData{ii}.timestep = 'unknown';
        end
    elseif isfield(sData{ii}, 'time')
        avgDiff = mode(diff(sort(sData{ii}.time)));
       
%         avgDiff = nanmean(diff(sort(sData{ii}.time)));
        if avgDiff < 1.5 && avgDiff > 0.5
            sData{ii}.timestep = 'daily';
        elseif avgDiff < 45 && avgDiff > 10
            sData{ii}.timestep = 'monthly';
        elseif avgDiff < 380 && avgDiff > 350
            sData{ii}.timestep = 'yearly';
        elseif isnan(avgDiff)
           sData{ii}.timestep = 'unknown';
        else
            error('methodDirect:unknownTimeStep',['The current time step is ' ...
                num2str(avgDiff) ', which has not been programmed for.']);
        end
    else
       sData{ii}.timestep = 'unknown';
    end
    warning('on', 'MATLAB:mode:EmptyInput')
    
    %Convert units:
%     sData{ii} = struct_2_standard_units(sData{ii}, varLd, sData{ii}.timestep);
    
    if ~isempty(unitsTable)
        nmAtt = ['att' varLd];
        unitsCurr = find_att(sData{ii}.(nmAtt), 'units');
        unitsUse  = find_att(unitsTable, varLd);
        
        if ~isequal(unitsCurr, unitsUse) && ~isempty(unitsUse)
            if regexpbl(unitsCurr, 'unknown')
                warning('dsLdFields:unknownUnits',['Units will not be switched because current '...
                    'units are unknown.']);
            else
                if strcmpi(unitsCurr, 'mm') && strcmpi(unitsUse, 'm') 
                    sData{ii}.(varLd) = sData{ii}.(varLd)/1000;
                    sData{ii}.(nmAtt) = set_att(sData{ii}.(nmAtt), 'units', unitsUse);
                elseif strcmpi(unitsCurr, 'm') && strcmpi(unitsUse, 'mm')
                    sData{ii}.(varLd) = sData{ii}.(varLd)*1000;
                    sData{ii}.(nmAtt) = set_att(sData{ii}.(nmAtt), 'units', unitsUse);
                else
                    warning('dsLdFields:unitsNotProgrammed', ['Units will ' ...
                        'not be switched because the combination has not ' ...
                        'been programmed for. Current units are ' unitsCurr  ...
                        ' and requested units are ' unitsUse '.']);
                end
            end
        end
    end
    
    if regexpbl(fldsLd{ii}, {'hist','sim'}, 'and')
        indSimHist(end+1) = ii;
    end
    if regexpbl(fldsLd{ii}, {'proj','sim'}, 'and')
        indSimProj(end+1) = ii;
    end
end


%If projection, but both projeciton and historic simulations present, stitch together
if blStitch == 1
%     warning('dsLdFlds:stitch','Ensure stitch edits work to match not only variables, but also sim4bc.')
    if ~isempty(indSimProj) && ~isempty(indSimHist) 
        if numel(indSimProj) == 1 && numel(indSimHist) == 1 
            varProjCurr = sData{indSimProj}.(varName);
            indUnd = regexpi(varProjCurr, '_');
            varProjCurr = varProjCurr(indUnd(1)+1:end);

            varHistCurr = sData{indSimHist}.(varName);
            indUnd = regexpi(varHistCurr, '_');
            varHistCurr = varHistCurr(indUnd(1)+1:end);

            if ~isequal(varProjCurr, varHistCurr)
                error('dsLdFlds:StichVar',['The projection variable is ' ...
                    varProjCurr ' and the historical variable is ' varHistCurr ...
                    '. They must be the same to stitch the simulations.']);
            end
            
            indHistVar = 1;
            
            if isequal(varProjCurr, varHistCurr)
                varAll = {varProjCurr};
            else
                error('dsLdFlds:varDiff',['The projection variable is ' ...
                    varProjCurr ' and the historic variable is ' varHistCurr ...
                    '. They should be the same.']);
            end
        elseif numel(indSimProj) == 2 && numel(indSimHist) == 2
            %Find variables and order in projection simulations:
            varAll = cell(0,1);
            for ii = 1 : numel(indSimProj)
                varProjCurr = sData{indSimProj(ii)}.(varName);
                indUnd = regexpi(varProjCurr, '_');
                varProjCurr = varProjCurr(indUnd(1)+1:end);
                if ~any(strcmpi(varAll, varProjCurr))
                   varAll{end+1} = varProjCurr; 
                end
            end

            if numel(varAll(:)) ~= 2
               error('dsLdFlds:notTwo',[num2str(numel(varAll(:))) ' unique variables were not found (Two are expected).']); 
            end

            %Confirm same variables present in historic and get order:
            indHistVar = nan(numel(varAll(:)), 1);
            for ii = 1 : numel(indSimHist)
                varHistCurr = sData{indSimHist(ii)}.(varName);
                indUnd = regexpi(varHistCurr, '_');
                varHistCurr = varHistCurr(indUnd(1)+1:end);
                if ~any(strcmpi(varAll, varHistCurr))
                   error('dsLdFlds:newVar', [varHistCurr ' variable is present in historic simulation but not projection simulation.']); 
                end

                indHistVar(ii) = find(strcmpi(varAll, varHistCurr) == 1);
            end
            clear ii
        else
            error('dsLdFlds:numberStichFlds',[num2str(numel(indSimProj)) ...
                ' projection simulations and ' num2str(numel(indSimHist)) ...
                ' historical simulation were identified. This combination '...
                'has not been programmed for.']);
        end   

        
        %Stitch variables:
        for ii = 1 : numel(indSimHist)
            indAsn = indSimHist(indHistVar(ii));
            indGet = indSimProj(ii);
            desc = sData{indAsn}.(varName);
            indUnd = regexpi(desc,'_');
            desc = ['projhist' desc(indUnd(1):end)];
            sData{indAsn}.(varName) = desc;
            
            indChange = find(indiceIn == indGet);
            if ~isempty(indChange)
               indiceIn(indChange) = indAsn; 
            end

            %Check that lat and lon are the same:
            if ~isequal(sData{indAsn}.(varLat), sData{indGet}.(varLat)) || ~isequal(sData{indAsn}.(varLon), sData{indGet}.(varLon))
                warning('dsLdFlds:diffCrd','The current fields cannot be stitched because their coordinates are different.');
%                 continue 
            end
            indUnd = regexpi(varAll{ii}, '_');
            if ~isempty(indUnd)
                varCurr = varAll{ii}(indUnd(end)+1:end);
            else
                varCurr = varAll{ii};
            end
                
            sData{indAsn}.(varCurr) = cat(1, sData{indAsn}.(varCurr), sData{indGet}.(varCurr));
            if isfield(sData{indAsn}, varTime)
                sData{indAsn}.(varTime) = cat(1, sData{indAsn}.(varTime), sData{indGet}.(varTime));
            end
            if isfield(sData{indAsn}, [varTime '_bnds'])
                szBnds = sData{indAsn}.([varTime '_bnds']);
                if szBnds(1) > szBnds(2)
                    indCat = 1;
                elseif szBnds(2) >= szBnds(1)
                    indCat = 2;
                end
                sData{indAsn}.([varTime '_bnds']) = cat(indCat, sData{indAsn}.([varTime '_bnds']), sData{indGet}.([varTime '_bnds']));
            end
            if isfield(sData{indAsn}, varDate)
                sData{indAsn}.(varDate) = cat(1, sData{indAsn}.(varDate), sData{indGet}.(varDate));
            end

            %Remove one datastructure and update indices
            sData(indGet) = [];

            indGetGt = find(indSimProj > indGet);
            if ~isempty(indGetGt)
                indSimProj(indGetGt) = indSimProj(indGetGt) - 1; 
            end
            indAsnGt = find(indSimHist > indGet);
            if ~isempty(indAsnGt)
                indSimHist(indAsnGt) = indSimHist(indAsnGt) - 1;
            end
            indDsGt = find(indiceIn > indGet);
            if ~isempty(indDsGt)
                indiceIn(indDsGt) = indiceIn(indDsGt) - 1;
            end
        end
        clear ii
    end
end

%Resample
if ~strcmpi(resample, 'none')
    nData = numel(sData(:));
            
    if regexpbl(resample, {'up','area','wgt'}, 'and')
        %Find all candidate datasets for upscaling and identify indice of lowest resolution:
        indResample = nan(0, 1);
        indLR = nan;
        deltaLR = 0.0001;
        for ii = 1 : nData
            if regexpbl(sData{ii}.(varName), {'ts'}) && regexpbl(sData{ii}.(varName), {'ref','sim'})
                indResample(end+1) = ii;
                
                deltaCurr = nanmean([nanmean(abs(diff(sData{ii}.(varLat)))), nanmean(abs(diff(sData{ii}.(varLon))))]);;
                if numel(indResample) == 1 || deltaCurr > deltaLR
                    indLR = ii;
                    deltaLR = deltaCurr;
                end
            end
        end
        clear ii
        
        %Upscale
        indResample = setdiff(indResample, indLR);
        for ii = 1 : numel(indResample)
            indUnd = regexpi(sData{indResample(ii)}.(varName), '_');
            var = sData{indResample(ii)}.(varName)(indUnd(end)+1:end);
            [sData{indLR}, sData{indResample(ii)}, nRe] = upscale_geodata_v3(sData{indLR}, sData{indResample(ii)}, var, resample, 1);
            if nRe ~= 2 && ~isnan(nRe)
               error('dsLdFld:wrongFldUpscale','The wrong field was upscaled.') 
            end
        end
        clear ii
    elseif regexpbl(resample, {'common','nearest'}, 'and')
        %Find common grid:
        latOut = [];
        lonOut = [];
        for ii = 1 : nData
            if regexpbl(sData{ii}.(varName), {'ts'}) && regexpbl(sData{ii}.(varName), {'ref','sim'})
                latOut = unique(union(latOut, sData{ii}.(varLat)));
                lonOut = unique(union(lonOut, sData{ii}.(varLon)));
            end
        end
        nLatOut = numel(latOut);
        nLonOut = numel(lonOut);
        [meshLonOut, meshLatOut] = meshgrid(lonOut, latOut);
        
        %Resample all to common grid:
        for ii = 1 : nData
            [meshLonIn, meshLatIn] = meshgrid(sData{ii}.(varLon), sData{ii}.(varLat));
            if regexpbl(sData{ii}.(varName), {'ts'}) && regexpbl(sData{ii}.(varName), {'ref','sim'})
                nTime = numel(sData{ii}.(varDate)(:,1));
                tempData = nan([nTime,nLatOut, nLonOut], 'single');
                
                indUnd = regexpi(sData{ii}.(varName), '_');
                var = sData{ii}.(varName)(indUnd(end)+1:end);
            
                for jj = 1 : nTime
                    tempData(jj,:,:) = griddata(meshLonIn, meshLatIn, double(squeeze(sData{ii}.(var)(jj,:,:))), meshLonOut, meshLatOut, 'nearest');
                end
                
                sData{ii}.(var) = tempData;
                sData{ii}.(varLat) = latOut;
                sData{ii}.(varLon) = lonOut;
            end
        end
        
    elseif regexpbl(resample, {'up','area','conserve','remap'}, 'and')    
        warning('dsLdFlds:unknownareaConserveRemap','area_conserve_remap has not yet been added. It should be consistent with the area_wgt option.')
    else
        error('dsLdFlds:unknownResampleMethod',['The resampling method ' ...
            resample ' has not been programmed for.']);
    end
end

if nargout > 1
    varargout{1} = indiceIn;
    if nargout > 2
        fldsTemp = cell(numel(sData(:)), 1);

        for ii = 1 : numel(sData(:))
            if isfield(sData{ii}, 'description')
                fldsTemp{ii} = sData{ii}.('description');
            end
        end
        varargout{2} = fldsTemp;
    end
end
