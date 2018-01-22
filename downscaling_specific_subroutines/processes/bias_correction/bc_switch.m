function [sDataOut, bcMethod] = bc_switch(dataIn, indSim2Bc, sDs)

if ~iscell(dataIn)
    error('bcSwitch:notCell','Input data are expcted to be a cell.');
end


%Choose bias correction time methodology based on time resolution of inputs
if isfield(dataIn{indSim2Bc}, 'timestep')
    if regexpbl(dataIn{indSim2Bc}.timestep, {'day', 'daily'})
        sDs.bctype = '15day_window';
    elseif regexpbl(dataIn{indSim2Bc}.timestep, {'month','yearly'})
        sDs.bctype = 'month';
    else
        error('methodDirect:tstepBc', ['The time step ' dataIn{indSim2Bc}.timestep ' has not been programmed for.']);
    end
else
    error('methodDirect:tstepUnknown', ['The time step for ' sDs.fldsDs  ' is unknown.']);
end

%Set type of bias correction based on variable:
if ~regexpbl(sDs.method, {'_mult', '_add'})
    bcScale = 'add';
    
%     if regexpbl(sDs.varDs, {'pr','rsds','rsdt'})
%         bcScale = 'add';
%     elseif regexpbl(sDs.varDs, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
%         bcScale = 'mult';
%     else
%         error('bcSwitch:unknownVar',['Variable ' sDs.varDs ' has not '...
%             'been programmed for. Determine which bias correction option to use.']);
%     end
    
    disp(['Autoselecting bias correction scaling method to be ' bcScale '.']);
else
    if regexpbl(sDs.method, '_mult')
        bcScale = 'mult';
    elseif regexpbl(sDs.method, '_add')
        bcScale = 'add';
    else
        bcScale = 'add';
        warning('bcSwitch:methodSet', ['Autoselecting bias correction scaling method to be ' bcScale '.']);
    end
end
   
sDs.bctype = [sDs.bctype '_' bcScale];


%Find which fields to use for bias correction (sim versus ref)
if isfield(sDs, 'bctype')
    type = sDs.bctype;
else
    warning('bcSwitch:noBcType', ['The field ' char(39) 'bctype' char(39) ' is not present.']);
    type = sDs.method;
end

bcMethod = 'unknown';

varNm = 'description';


if ~regexpbl(sDs.method,'eQM') && regexpbl(sDs.method, '2var')
    error('bcSwitch:bcMethdNotProgrammed', ['The downscaling method '...
        'indicates two variables to be used in bias correction but the '...
        'method is not QM. Currently the proxy variable method can only be used with QM.']);
end


%Assign input data to bias correction inputs
for ii = 1 : numel(dataIn(:))
    nmDataCurr = dataIn{ii}.(varNm);
    
    %Find simulation baseline dataset:
    if regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, {'sim','hist', '4bc'}, 'and')
        sSimBase = dataIn{ii};
        indUnd = regexpi(nmDataCurr, '_');
        varBase = nmDataCurr(indUnd(end)+1:end);
    elseif regexpbl(nmDataCurr, {sDs.varDs, 'sim','hist'}, 'and') && ~regexpbl(nmDataCurr, '4bc')
        sSimBase = dataIn{ii};
        indUnd = regexpi(nmDataCurr, '_');
        varBase = nmDataCurr(indUnd(end)+1:end);
    end
    
    %Find simulation dataset to bias correct:
    if ~regexpbl(nmDataCurr, 'ref') && regexpbl(nmDataCurr, sDs.varDs) && ((regexpbl(sDs.period, 'proj') && regexpbl(nmDataCurr, 'proj')) || (regexpbl(sDs.period, 'hist') && regexpbl(nmDataCurr, 'hist')))
        sSim2Bc = dataIn{ii};
        indUnd = regexpi(nmDataCurr, '_');
        var2Bc = nmDataCurr(indUnd(end)+1:end);
    end
    
    %Find reference dataset to be used in bias correction:
    if regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, {'ref','hist'}, 'and')
        sRefBase = dataIn{ii};
        indUnd = regexpi(nmDataCurr, '_');
        varRef = nmDataCurr(indUnd(end)+1:end);
    elseif ~regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, {sDs.varDs, 'ref','hist'}, 'and')
        sRefBase = dataIn{ii};
        indUnd = regexpi(nmDataCurr, '_');
        varRef = nmDataCurr(indUnd(end)+1:end);
    end
    
    
%     if regexpbl(nmDataCurr, {sDs.varDs, 'sim'}, 'and')
%         if regexpbl(sDs.period, 'proj')
%             if regexpbl(nmDataCurr, 'proj')
%                 sSim2Bc = dataIn{ii};
%                 indUnd = regexpi(nmDataCurr, '_');
%                 var2Bc = nmDataCurr(indUnd(end)+1:end);
%             end
%             
%             if regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, {'hist', '4bc'}, 'and')
%                 sSimBase = dataIn{ii};
%                 indUnd = regexpi(nmDataCurr, '_');
%                 varBase = nmDataCurr(indUnd(end)+1:end);
%             elseif ~regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, 'hist')
%                 sSimBase = dataIn{ii};
%                 indUnd = regexpi(nmDataCurr, '_');
%                 varBase = nmDataCurr(indUnd(end)+1:end);
%             end
%         elseif regexpbl(sDs.period, 'hist')
%             sSim2Bc = dataIn{ii};
%             indUnd = regexpi(nmDataCurr, '_');
%             var2Bc = nmDataCurr(indUnd(end)+1:end);
%             if regexpbl(sDs.method, '2var') && regexpbl(nmDataCurr, '4bc')
%                 sSimBase = dataIn{ii};
%                 indUnd = regexpi(nmDataCurr, '_');
%                 varBase = nmDataCurr(indUnd(end)+1:end);
%             elseif ~regexpbl(sDs.method, '2var')
%                 sSimBase = dataIn{ii};
%                 indUnd = regexpi(nmDataCurr, '_');
%                 varBase = nmDataCurr(indUnd(end)+1:end);
%             end
%         else
%             error('bcSwitch:projPeriod', [sDs.period ' is not a recognized downscaling period.'])
%         end
%     elseif regexpbl(nmDataCurr, {sDs.varDs, 'ref','hist'}, 'and')
%         sRefBase = dataIn{ii};
%         indUnd = regexpi(nmDataCurr, '_');
%         varBase = nmDataCurr(indUnd(end)+1:end);
%     end
    
    if sum(~isnan(sDs.indRef)) == 2
        if ~regexpbl(nmDataCurr, sDs.varDs) %Current variable different from variable being downscaled
            indUnd = regexpi(nmDataCurr, '_');
            if ~isempty(indUnd)
                varV = nmDataCurr(indUnd(end)+1:end);
            else
                varV = nmDataCurr;
            end
            if regexpbl(nmDataCurr, 'sim')
                if regexpbl(sDs.period, 'proj')
                    if regexpbl(nmDataCurr, 'proj')
                        sSim2BcV = dataIn{ii};
                    end
                    if regexpbl(nmDataCurr, 'hist')
                        sSimBaseV = dataIn{ii};
                    end
                elseif regexpbl(sDs.period, 'hist')
                    sSim2BcV = dataIn{ii};
                    sSimBaseV = dataIn{ii};
                else
                    error('bcSwitch:projPeriod', [sDs.period ' is not a recognized downscaling period.'])
                end
            elseif regexpbl(nmDataCurr, {'ref','hist'}, 'and')
                sRefBaseV = dataIn{ii}; 
            end
        end
    end
end


%Check that all inputs have same time resolution:
if ~isequal(sSim2Bc.timestep, sSimBase.timestep) || ~isequal(sSim2Bc.timestep, sRefBase.timestep)
    error('bcSwitch:diffTimeStep','The inputs to bias correct have difference time resolutions. Therefore, bias correction cannot proceed.')
end

%ENSURE VARIABLES TO USE IN BIAS CORRECTION MATCH THE DOWNSCALING VARIABLE
if ~isequal(varBase, sDs.varDs)
    temp = sSimBase.(varBase);
    sSimBase.(sDs.varDs) = temp;
    sSimBase = rmfield(sSimBase, varBase);
end
if ~isequal(varRef, sDs.varDs)
    temp = sRefBase.(varRef);
    sRefBase.(sDs.varDs) =  temp;
    sRefBase = rmfield(sRefBase, varRef);
end

if ~isequal(var2Bc, sDs.varDs)
    error('bcSwitch:diffVar',['The variable chosen for bias correction is ' ...
        var2Bc ' but the downscaling variable is ' sDs.varDs '. This is not expected.'])
end


%IMPLEMENT BIAS CORRECTION METHODS
blJnt = 0;
nPtsPlt = 2;
[pathCdfRt, ~, ~] = fileparts(sDs.pathOut);
pathCDF = fullfile(pathCdfRt,'figures','cdf');
if regexpbl(sDs.method, {'eQM','q2q'}) %Use empirical quantile mapping
%     sSim2Bc
%     sSimBase
%     sRefBase

    if regexpbl(sDs.method, 'eQM')
        bcMethod = 'eQM';
    elseif regexpbl(sDs.method, 'q2q')
        bcMethod = 'q2q';
    else
        error('bcSwitch:bcType1D',['The method, ' sDs.method ', does not contain a known bias correction method']);
    end
    
    type = [bcMethod, '_', type];
    
    sDataOut = bc1d_geodata(sSim2Bc, sSimBase, sRefBase, sDs.varDs, sDs.yrsBase, type);
    
elseif regexpbl(sDs.method,'eJBC') %Use empirical joint bias correction
    if any(isnan(sDs.indRef)) || numel(sDs.indRef) ~= 2
       error('bcSwitch:wrongNumbRefJbc',['JBC bias correction requires two '...
           'reference dataset, but ' num2str(numel(sDs.indRef)) ' are present.']); 
    end
    
    if regexpbl(sDs.method, 'eJBC')
        bcMethod = 'eJBC';
    else
        error('bcSwitch:bcType2D',['The method, ' sDs.method ', does not contain a known bias correction method']);
    end
        
    type = [bcMethod, '_', type];
    
    %Bias-correct time-series being downscaled:
    [sDataOut, sDataOutV] = bc2d_geodata(sSim2Bc, sSim2BcV, sSimBase, sSimBaseV, sRefBase, sRefBaseV, sDs.varDs, varV, sDs.yrsBase, type);
    
    if regexpbl(varV, {'pr','rsds','rsdt'})
        sDataOutV.(varV)(sDataOutV.(varV) < 0) = 0;
    elseif ~regexpbl(varV, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
        error('bcSwitch:unknownVar',['Variable ' sDs.varDs ' has not '...
            'been programmed for. Determine which bias correction option to use.']);
    end

    blJnt = 1;
end


%Ensure specific variables are not negative
if regexpbl(sDs.varDs, {'pr','rsds','rsdt'})
    sDataOut.(sDs.varDs)(sDataOut.(sDs.varDs) < 0) = 0;
elseif ~regexpbl(sDs.varDs, {'tas','tmp','tav','tmax', 'tmx', 'tmn', 'tmin','temp'})
    error('bcSwitch:unknownVarNeg',['Variable ' sDs.varDs ' has not '...
        'been programmed for. Determine the domain for this variable.']);
end


%%Make plot of joint histograms if requested:
if isfield(sDs, 'wrtOut') && regexpbl(sDs.wrtOut,'cdf')
    if blJnt == 1
        print_hist2d(pathCDF, sDataOut, sDataOutV, sSim2Bc, sSim2BcV, sSimBase, sSimBaseV, sRefBase, sRefBaseV, sDs.varDs, varV, sDs.yrsDs, sDs.yrsBase, nPtsPlt);
    else
        print_CDF(pathCDF, sDataOut, sSim2Bc, sSimBase, sRefBase, sDs.varDs, sDs.yrsDs, sDs.yrsBase, nPtsPlt);
    end
end


