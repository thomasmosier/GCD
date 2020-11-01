function method_extract(sPath, sDs)


%Set select names of variables for output structure array:
varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';
varData = sDs.varDs;


%Make figures
blWrtFig = 1;

%Path for storing output figures:
if blWrtFig == 1
    foldPlot = fullfile(sPath.output, 'figures');
    if ~exist(foldPlot, 'dir')
        mkdir(foldPlot)
    end
    
    
    if regexpbl(varData, {'tas','tmp','tave','tavg', 'tmn', 'tmin', 'tmx', 'tmax'})
        typePlot = 'average';
    elseif regexpbl(varData, {'pre','rsds'}) || strcmpi(varData, 'pr')
        typePlot = 'total';
    else
        error('methodExtract:unknownVar',['Variable ' varData ' has not been programmed for.'])
    end
end
warning('off', 'MATLAB:LargeImage');

varNm = 'Name';

%Load station coordinate file:
if numel(sDs.fldsTInv(:)) == 1
    [~, ~, ext] = fileparts(sPath.(sDs.fldsTInv{1}));
    if regexpbl(ext,'xls')
        [~, ~, raw] = xlsread(sPath.(sDs.fldsTInv{1}));
        crdFlds = raw(1,:);
        raw = raw(2:end,:);
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
        
        sCrd = struct;
        if any(strcmpi(crdFlds, varLat)) && any(strcmpi(crdFlds, varLon)) && any(strcmpi(crdFlds, varNm))
            indLat = find(strcmpi(crdFlds, varLat) == 1);
            indLon = find(strcmpi(crdFlds, varLon) == 1);
            indNm  = find(strcmpi(crdFlds,  varNm) == 1);
            
            sCrd(:).(varLat) = cell2mat(raw(:,indLat));
            sCrd(:).(varLon) = cell2mat(raw(:,indLon));
            sCrd(:).( varNm) = raw(:, indNm);
            
            crdFlds = {varNm, varLat, varLon};
            nmStn = varNm;
        else
            for ii = 1 : numel(crdFlds)
    %             sCrd = setfield(sCrd, crdFlds{ii}, raw(:,ii));
                if ischar(raw(1,ii)) || iscell(raw(1,ii))
                    if iscell(raw{1,ii})
                        colTemp = cell(numel(raw(:,ii)), 1);
                        for jj = 1 : numel(raw(:,ii))
                            colTemp{jj} = char(raw{jj,ii});
                        end
                        sCrd(:).(crdFlds{ii}) = colTemp;
                    elseif isnumeric(raw{1,ii})
                        sCrd(:).(crdFlds{ii}) = cell2mat(raw(:,ii));
                    else
                        sCrd(:).(crdFlds{ii}) = raw(:,ii);
                    end
                elseif isnumeric(raw(1,ii))
                    sCrd(:).(crdFlds{ii}) = cell2mat(raw(:,ii));
                else
                    error('methodExtract:stnFldFormat',['The current column '...
                        'is class ' class(raw(1,ii)) ', which has not been programmed for.']);
                end
            end
            clear ii raw
        end
        
        %Find lat and lon fields
        nmStnLat = blanks(0);
        nmStnLon = blanks(0);
        nmStn = blanks(0);
        for ii = 1 : numel(crdFlds)
            if regexpbl(crdFlds{ii}, 'lat')
                nmStnLat = crdFlds{ii};
            elseif regexpbl(crdFlds{ii}, 'lon')
                nmStnLon = crdFlds{ii};
            elseif regexpbl(crdFlds{ii}, 'name')
                nmStn = crdFlds{ii};
            end
        end
        clear ii
        
        %Remove spaces and cahracters:
        strrep(nmStn, ' ', '');
        strrep(nmStn, '-', '');
        strrep(nmStn, '_', '');
        

        %Find station coordinate bounds (used for loading any reference data):
        sDs.lonDs = [min(sCrd.(nmStnLon)), max(sCrd.(nmStnLon))]; %[West, East] (in decimal degrees; -180 to 180)
        sDs.latDs = [min(sCrd.(nmStnLat)), max(sCrd.(nmStnLat))]; %[North, South] (in decimal degrees; 90 to -90)
        
        %Number of stations:
        nStn = numel(sCrd.(crdFlds{1})(:));
    else
        error('methodExtract:stnCrdFormat',['The current station ' ...
            'coordinate file is of type ' ext '. Excel format is required.']);
    end
else
    error('methodExtract:stnCrdNotPresent',[num2str(numel(sDs.fldsTInv(:)))  ...
        ' time invariant files present. Exactly one is expected.']);
end


%Load station data file:
if regexpbl(sDs.method, 'stn')
    if numel(sDs.fldsTVar(:)) >= 1
        indStnData = nan(0,1);
        for ii = 1 : numel(sDs.fldsTVar(:))
            if regexpbl(sDs.fldsTVar{ii}, 'rec')
                indStnData(end+1) = ii; 
            end
        end
        clear ii

        if numel(indStnData) ~= 1
            error('methodExtract:stnDataNotPresent', [num2str(numel(indStnData))  ...
                ' time varying files present. Exactly one is expected.']);
        end

        pathStnData = sPath.(sDs.fldsTVar{indStnData});
        [~, ~, raw] = xlsread(pathStnData);
        raw = raw(2:end,:);
        raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

        % Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
        raw(R) = {NaN}; % Replace non-numeric cells

        raw = reshape([raw{:}],size(raw));

        %Find header indices and reformat dates/data
        [~, stnDataNm] = xlsread(pathStnData);
        stnDataNm = stnDataNm(1,:);
        
        %Remove 'date' if present:
        indStnDate = find(strcmpi(stnDataNm, 'date') == 1);
        if ~isempty(indStnDate)
            stnDataNm(indStnDate) = [];
            raw(:, indStnDate) = [];
        end
        
        %Find date components:
        indStnYr = find(strcmpi(stnDataNm, 'year') == 1);
        indStnMnth = find(strcmpi(stnDataNm, 'month') == 1);
        indStnDy = find(strcmpi(stnDataNm, 'day') == 1);
        indStnHdr = unique([indStnYr, indStnMnth, indStnDy]);
        if numel(indStnHdr) ~= 3
            error('methodExtract:unknownHdr', [num2str(numel(indStnHdr)) ...
                ' header elements found, but three expected in ' pathStnData]);
        else
            stnDataNm(indStnHdr) = [];
        end

        stnDate = raw(:, [indStnYr, indStnMnth, indStnDy]);
        stnData = raw(:, numel(indStnHdr) + 1 : end);
    else
        error('methodExtract:tVarDataNotPresent',[num2str(numel(sDs.fldsTVar(:)))  ...
            ' time varying files present. At least one is expected.']);
    end
    
    %Determine which time varying fields are stations vs. grids:
    fldStn = cell(0,1);
    for ii = numel(sDs.fldsTVar) : -1 :1
        if regexpbl(sDs.fldsTVar{ii}, {'rec','stn'})
            fldStn{end+1} = sDs.fldsTVar{ii};
            sDs.fldsTVar(ii) = [];
        end
    end

    if numel(fldStn(:)) ~= 1
       error('methodExtract:noStnFld','No station field found in cell array of time-varying fields.'); 
    end
else
    stnDataNm = sCrd.(nmStn)(:)';
    stnDate = nan(1,3);
    stnData = nan(1,nStn);
end


%Check that there are same number/order of coordinates and stations; assign
%to single structure array
sCrd.(varData) = cell(nStn,1);
sCrd.(varDate) = cell(nStn,1);

if numel(stnDataNm(:)) ~= numel(sCrd.(crdFlds{1})(:))
    error('methodExtract:stnNumber', [num2str(numel(stnData(1,:)))  ...
        ' stations found in data file and ' num2str(numel(sCrd(:).(crdFlds{1})(:))) ...
        ' found in coordinate file. These numbers are expected to be the same.']);
else
    indOrd = (1:nStn);
    for ii = nStn : -1 : 1
        %Find correct index
        indMatchCurr = find(strcmpi(sCrd.(nmStn){ii}, stnDataNm(indOrd)) == 1);
        if isempty(indMatchCurr)
            error('methodExtract:stnNotPresent',[sCrd.(nmStn){ii} ...
                ' is present in the coordinate file but appears to be' ...
                ' missing from input station records.']);
        elseif indMatchCurr ~= ii %Switch order
            indTemp = indOrd(ii);
            indOrd(ii) = indMatchCurr;
            indOrd(indMatchCurr) = indTemp;
        end
        
        %Assign to structure
        sCrd.(varData){ii} = stnData(:, indOrd(ii));
        sCrd.(varDate){ii} = stnDate;
        
        %Keep only non-nan dates/data
        indKp = find(~isnan(sCrd.(varData){ii}));
        sCrd.(varData){ii} = sCrd.(varData){ii}(indKp);
        sCrd.(varDate){ii} = sCrd.(varDate){ii}(indKp,:);
    end
    clear ii
    
    %Check that swapped names match coordinate structure
    if ~isequal(sCrd.(nmStn), stnDataNm(indOrd)')
        error('methodExtract:stnOrder', ['The station order does not '...
            'match between the data and coordinate files.']);
    end
end


%Set downscaling indice:
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

%Get original write format (needed because it is switched to csv for
%points:
orgWrtTyp = sDs.wrtFmt;

warning('off', 'ncCal:unknownCal'); 
warning('off', 'eQMgeodata:unequalYrs');


tic;
%Loop over months (if timestep is monthly)
for ii = 1 : sDs.nLp
    if regexpbl(sDs.timestep, {'month', 'mnth'})
        mnthCurr = sDs.mnthsDs(ii);
        mnthDisp = mnth_str(ii);
    elseif regexpbl(sDs.timestep, {'day', 'daily'})
        mnthCurr = nan;
        mnthDisp = 'annual';
    end

    %%LOAD NON-STATION INPUTS:
        indDsOrg = sDs.indDs;
    [sTVar, indNew] = ds_ld_fields(sPath, sDs.fldsTVar, sDs.lonDs, sDs.latDs, sDs.yrsLd, mnthCurr, 'units', sDs.units, 'stitch', 1, 'indice', [sDs.indDs(:); sDs.indRef(:)], 'resample', sDs.resample, 'frame', 2);
    nTVar = numel(sTVar);
    %Update 'indDsIn' based on change in indDs
        sDs.indDs = indNew(1);
        sDs.indRef = indNew(2:end);
        indDelta = (sDs.indDs - indDsOrg);
        indDsIn = indDsIn + indDelta;
    
    %loop over stations
    for kk = 1 : nStn
        %Reset downscaling inde for each station:
        sDs.indDs = indDsIn;
        
        %Create structure array for current station:
        sStnCurr = struct;
            sStnCurr.(varLat) = sCrd.(nmStnLat)(kk);
            sStnCurr.(varLon) = sCrd.(nmStnLon)(kk);
            sStnCurr.(varDate) = sCrd.(varDate){kk};
            sStnCurr.(varData) = sCrd.(varData){kk};
        
        %Write station records to file:
        if ~isempty(sStnCurr.(varData))
                sDs.wrtFmt = 'csv';
            ds_wrt_outputs(sStnCurr, sCrd.(nmStn){kk}, sDs, sPath, 'fold', 'stn_records');
                sDs.wrtFmt = orgWrtTyp;
        end
        
        %Create structure array for all time-varying inputs at current
        %station point
        sTVarPt = cell(nTVar, 1);
        %Find indices in gridded data for current station
        indTVarLat = nan(nTVar,1);
        indTVarLon = indTVarLat;
        for ll = 1 : nTVar
            [~, indTVarLat(ll)] = min(abs(sTVar{ll}.(varLat) - sStnCurr.(varLat)));
            [~, indTVarLon(ll)] = min(abs(sTVar{ll}.(varLon) - sStnCurr.(varLon)));
            
            sTVarPt{ll} = sTVar{ll};
                sTVarPt{ll}.(varLat) = sTVar{ll}.(varLat)(indTVarLat(ll));
                sTVarPt{ll}.(varLon) = sTVar{ll}.(varLon)(indTVarLon(ll));
                sTVarPt{ll}.(sDs.varDs) = sTVar{ll}.(sDs.varDs)(:, indTVarLat(ll), indTVarLon(ll));
                
            %Write inputs at each location:
            %Get base file name:
            nmBaseOut = ds_file_nm(sPath, sDs, ll);

            %Write station specific outputs:
                sDs.wrtFmt = 'csv';
            ds_wrt_outputs(sTVarPt{ll}, sCrd.(nmStn){kk}, sDs, sPath, 'fold', nmBaseOut);
                sDs.wrtFmt = orgWrtTyp;
        end
        clear ll
        
        %Get base file name:
        nmBaseOut = ds_file_nm(sPath, sDs, indDsIn);
        
        %If 'stnAvg' bias correction is programmed for, note and remove
        %from method string (occurs after any other bias correction method:
        blStnAvgBc = 0;
        blStnQmBc = 0;
        methodFull = sDs.method;
        if regexpbl(sDs.method, 'stnAvg')
           blStnAvgBc = 1; 
           sDs.method = strrep(strrep(sDs.method, 'stnAvg', ''), '__', '_');
           if isequal(sDs.method(1), '_')
               sDs.method = sDs.method(2:end);
           end
        elseif regexpbl(sDs.method, 'stnQM')
           blStnQmBc = 1; 
           sDs.method = strrep(strrep(sDs.method, 'stnQM', ''), '__', '_');
           if isequal(sDs.method(1), '_')
               sDs.method = sDs.method(2:end);
           end
        end
        
        %Bias Correct gridded data (before removing station bias):
        if ~isempty(sDs.indRef) && ~isnan(sDs.indRef) 
            [sTVarPt{end+1}, bcMethod] = bc_switch(sTVarPt, indDsIn, sDs);
            
            sDs.fldsTVar{numel(sTVarPt(:))} = [sDs.fldsTVar{indDsIn} '_bc-' bcMethod];
            sDs.indDs = numel(sDs.fldsTVar);
            
                sDs.wrtFmt = 'csv';
            ds_wrt_outputs(sTVarPt{sDs.indDs}, [sCrd.(nmStn){kk} '-bc-' bcMethod], sDs, sPath, 'fold', nmBaseOut, 'yrs', sDs.yrsDs);
                sDs.wrtFmt = orgWrtTyp;
        else 
            sDs.indDs = indDsIn;
        end
        
        sDs.method = methodFull;
        
        %Implement station average bias removal
        if blStnAvgBc == 1 || blStnQmBc == 1
            sTVarPtStnBc = sTVarPt;
            
            if blStnAvgBc == 1
                for ll = 1 : numel(sDs.fldsTVar)
                    if ll ~= sDs.indDs
                        warning('off', 'stnAvgBc:noStnRec');
                    end
                    sTVarPtStnBc{ll} = stn_avg_bc(sStnCurr, sTVarPt{ll}, sSim2Bc, sDs.varDs);
                    if ll ~= sDs.indDs
                        warning('on', 'stnAvgBc:noStnRec');
                    end
                end
                clear ll
            elseif blStnQmBc == 1
                %Find which fields to use for bias correction (sim versus ref)
                sTVarPtStnBc{1} = bc_switch_stn(sStnCurr, sTVarPt, indDsIn, sDs);
                %sTVarPtStnBc{ll} = stn_eQM_bc(sStnCurr, sTVarPt{ll}, sSim2Bc, sDs.varDs, sDs.timestep);
            end
              
            
            %Write station specific outputs:
                sDs.wrtFmt = 'csv';
            %(Bias corrected with station only)
            if blStnAvgBc == 1
                nmBcStnOnly = '-bc-avg-stn-only';
                if ~exist('bcMethod', 'var')
                   bcMethod = sDs.method; 
                end
                nmBcComb = ['-bc-' bcMethod '-global-and-avg-stn'];
            elseif blStnQmBc == 1
                nmBcStnOnly = '-bc-eQM-stn-only';
                if ~exist('bcMethod', 'var')
                   bcMethod = sDs.method; 
                end
                nmBcComb = ['-bc-' bcMethod '-global-and-eQM-stn'];
            end
            
            %Bias corrected with station data only
            ds_wrt_outputs(sTVarPtStnBc{indDsIn}, [sCrd.(nmStn){kk} nmBcStnOnly], sDs, sPath, 'fold', nmBaseOut, 'yrs', sDs.yrsDs);    
            %Bias corrected with both the global method and station records
            ds_wrt_outputs(sTVarPtStnBc{sDs.indDs}, [sCrd.(nmStn){kk} nmBcComb], sDs, sPath, 'fold', nmBaseOut, 'yrs', sDs.yrsDs);
                sDs.wrtFmt = orgWrtTyp;
        end
        
        %Plot annual and monthly time-series
        if blWrtFig == 1
            dataPlot = cell(numel(sTVarPt(:)), 1);
            datesPlot = dataPlot;
            for ll = 1 : numel(sTVarPt(:))
                dataPlot{ll} = sTVarPt{ll}.(varData);
                datesPlot{ll} = sTVarPt{ll}.(varDate);
            end
            
            %Get units:
            try
                unitsPlot = find_att(sTVarPt{sDs.indDs}.(['att' sDs.varDs]), 'units');
            catch
                unitsPlot = nan;
            end
            
            pathAnnPlot  = fullfile(foldPlot, [sCrd.(nmStn){kk}  '_annual_ts']);
            pathMnthPlot = fullfile(foldPlot, [sCrd.(nmStn){kk} '_monthly_ts']);
            
            plot_ts( pathAnnPlot, dataPlot, datesPlot, sDs.fldsTVar, sDs.yrsDs, 'annual', typePlot, sDs.varDs, unitsPlot);
            plot_ts(pathMnthPlot, dataPlot, datesPlot, sDs.fldsTVar, sDs.yrsDs,  'month', typePlot, sDs.varDs, unitsPlot);
        end
        
        %Gather all station outputs (write to single file outside of station loop)
        %Initialize output cell array
        if kk == 1
            sPtAllOut = cell(numel(sTVarPt), 1);
            
            if blStnAvgBc == 1 || blStnQmBc == 1
                sPtAllOutStnBc = cell(numel(sTVarPt), 1);
            end
        end
        for ll = 1 : numel(sTVarPt)
            if kk == 1
                sPtAllOut{ll} = struct;
                sPtAllOut{ll}.(varDate) = sTVarPt{ll}.(varDate);
                if blStnAvgBc == 1 || blStnQmBc == 1
                    sPtAllOutStnBc{ll} = struct;
                    sPtAllOutStnBc{ll}.(varDate) = sTVarPt{ll}.(varDate);
                end
            end
            sPtAllOut{ll}.(sCrd.(nmStn){kk}) = sTVarPt{ll}.(sDs.varDs);
            
            if blStnAvgBc == 1 || blStnQmBc == 1
                sPtAllOutStnBc{ll}.(sCrd.(nmStn){kk}) = sTVarPtStnBc{ll}.(sDs.varDs);
            end
        end
    end %End loop over stations
    clear kk

    
    %Write gathered outputs (All inputs, including those that have been bias corrected):
    for ll = 1 : numel(sPtAllOut)
        if regexpbl(sDs.fldsTVar{ll}, 'hist')
            yrsWrt = sDs.yrsBase;
        elseif regexpbl(sDs.fldsTVar{ll}, 'proj')
            yrsWrt = sDs.yrsDs;
        else
            yrsWrt = nan(1,2);
        end
        
        nmBaseOutCurr = ds_file_nm(sPath, sDs, ll);
        pathDataCurr = fullfile(sPath.output, 'all_pts', [nmBaseOutCurr '_all-pts']);
        
        wrt_csv_ts(pathDataCurr, sPtAllOut{ll}, sPtAllOut{ll}.(varDate), sCrd.(nmStn), 'yrs', yrsWrt);
    end
    clear ll
    
    
    %Write gathered outputs bias corrected with station)
    if blStnAvgBc == 1 || blStnQmBc == 1
        for ll = 1 : numel(sPtAllOut)
            nmBaseOutCurr = ds_file_nm(sPath, sDs, ll);
            if ll == sDs.indDs %(Bias corrected with previous method and station)
                nmFileOut = [nmBaseOutCurr '_all-pts' nmBcComb];
            else %(Bias corrected with station only)
                nmFileOut = [nmBaseOutCurr '_all-pts' nmBcStnOnly];
            end
            
            if regexpbl(sDs.fldsTVar{ll}, 'hist')
                yrsWrt = sDs.yrsBase;
            elseif regexpbl(sDs.fldsTVar{ll}, 'proj')
                yrsWrt = sDs.yrsDs;
            else
                yrsWrt = nan(1,2);
            end
            
            pathDataCurr = fullfile(sPath.output, 'all_pts', nmFileOut);
            wrt_csv_ts(pathDataCurr, sPtAllOutStnBc{ll}, sPtAllOutStnBc{ll}.(varDate), sCrd.(nmStn), 'yrs', yrsWrt);
        end
        clear ll
    end
    
%         fldsAll = sDs.fldsTVar;
%         sDs.fldsTVar = sCrd.(nmStn);
%         sDs.wrtFmt = 'csv';
%     ds_wrt_inputs(sPtAllOut, sDs, sPath);
%     if blStnAvgBc == 1
%         ds_wrt_inputs(sPtAllOutStnBc, sDs, sPath);
%     end
%         sDs.fldsTVar = fldsAll;
%         sDs.wrtFmt = orgWrtTyp;
    
    %Write inputs:
    if ~isempty(sDs.wrtOut)
            fldsAll = sDs.fldsTVar;
            sDs.fldsTVar = sDs.fldsTVar(1:numel(sTVar(:)));
        ds_wrt_inputs(sTVar, sDs, sPath);
            sDs.fldsTVar = fldsAll;
    end

    %%DISPLAY DURATION OF TIME FOR-LOOP HAS BEEN RUNNING
    deltatLoop = toc;
    perCmplt = 100*ii / numel(sDs.nLp);
    disp(['The ' sDs.method ' method has finished processing ' ...
        sDs.varDs ' data for ' mnthDisp ' (' ...
        num2str(round(100*perCmplt)/100) '% complete; elapsed time = ' ...
        num2str( round2(deltatLoop /60, 1)) ' minutes).' char(10)]);
end
clear ii

warning('on', 'ncCal:unknownCal');
warning('on', 'eQMgeodata:unequalYrs');
warning('on', 'MATLAB:LargeImage');