function sSimOut = stn_eQM_bc(sStnRec, sSimHist, sSim2Bc, varDs, type)



varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';

typeDs = 'eQM';



if ~regexpbl(typeDs, {'_mult', '_add'})
    typeDs = [typeDs '_add'];
    
%     if regexpbl(sDs.varDs, {'pr','rsds','rsdt'})
%         bcScale = 'add';
%     elseif regexpbl(sDs.varDs, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
%         bcScale = 'mult';
%     else
%         error('bcSwitch:unknownVar',['Variable ' sDs.varDs ' has not '...
%             'been programmed for. Determine which bias correction option to use.']);
%     end
    
    disp(['Autoselecting bias correction scaling method to be ' typeDs '.']);
else
    if regexpbl(typeDs, '_mult')
        typeDs = [typeDs '_mult'];
    elseif regexpbl(typeDs, '_add')
        typeDs = [typeDs '_add'];
    else
        typeDs = [typeDs '_add'];
        warning('bcSwitch:methodSet', ['Autoselecting bias correction scaling method to be ' typeDs '.']);
    end
end


%Find time-series elements to compare during bias correction:
%For historic model:
[indSimH, ~] = ismember(floor(sSimHist.(varDate)(:,1:end)), floor(sStnRec.(varDate)(:,1:end)), 'rows');
[indRefH, ~] = ismember(floor(sStnRec.(varDate)(:,1:end)), floor(sSimHist.(varDate)(:,1:end)),  'rows');

indSimH = find(indSimH ~= 0);
indRefH = find(indRefH ~= 0);

if numel(indSimH) ~= numel(indSimH)
    error('stnAvgBc:diffNumber',['Two different numbers of time-series ' ...
        'elements were selected. This is unexpected.']);
end


sSimOut = sSimHist;
if numel(indSimH) == 0 || numel(indRefH) == 0 
    warning('stnAvgBc:noStnRec',['No station records were found for times '...
        'overlapping with the simulation. Therefore the output will not be edited.']);
    return
else
    %Find indices for each output grid cell:
    [~, indLatSim] = min(abs(sStnRec.(varLat) - sSimHist.(varLat)));
    [~, indLonSim] = min(abs(sStnRec.(varLon) - sSimHist.(varLon)));

    if regexpbl(type, 'month')
        sSimOut.(varDs)(:, indLatSim, indLonSim) = bc1d_month(...
            sStnRec.(varDs)(indRefH), sStnRec.(varDate)(indRefH,:), ...
            sSimHist.(varDs)(indSimH,indLatSim,indLonSim), sSimHist.(varDate)(indSimH,:), ...
            sSim2Bc.(varDs)(:,indLatSim,indLonSim), sSim2Bc.(varDate)(:,:), ...
            typeDs);
    elseif regexpbl(type, {'day','window'})
        %Find number of days (i.e. +/- current)
        nDyWin = regexpi(type,'[\d]{1,3}day','match');
%             nDyWin = regexp(type,'-?\d+\.?\d*|-?\d*\.?\d+','match');
        if isempty(nDyWin) 
            nDyWin = 15;
        elseif iscell(nDyWin)
            if numel(nDyWin) > 1
                error('bc1d:multDayWindows',[num2str(numel(nDyWin)) ' day window indicators were found.'])
            else
                nDyWin = str2double(char(regexpi(nDyWin{1},'[\d]{1,3}','match')));
            end
        end
        
        sSimOut.(varDs)(:, indLatSim, indLonSim) = bc1d_day(...
            sStnRec.(varDs)(indRefH), sStnRec.(varDate)(indRefH,:), ...
            sSimHist.(varDs)(indSimH,indLatSim,indLonSim), sSimHist.(varDate)(indSimH,:), ...
            sSim2Bc.(varDs)(:,indLatSim,indLonSim), sSim2Bc.(varDate)(:,:), ...
            nDyWin, typeDs);
    else
        error('eQmGeodata:unknownType',['eQM type ' type ' has not been programmed for.']);
    end
end

if strcmpi(varDs, 'pre') || strcmpi(varDs, 'pr') 
    sSimOut.(varDs)(sSimOut.(varDs) < 0) = 0;
end