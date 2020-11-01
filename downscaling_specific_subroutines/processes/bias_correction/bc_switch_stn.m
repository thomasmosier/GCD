function sDataOut = bc_switch_stn(sStnCurr, dataIn, indSim2Bc, sDs)


if ~iscell(dataIn)
    error('bcSwitch:notCell','Input data are expcted to be a cell.');
end

varNm = 'description';


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


%Assign input data to bias correction inputs
for ii = 1 : numel(dataIn(:))
    nmDataCurr = dataIn{ii}.(varNm);
    
    %Find simulation baseline dataset:
    if regexpbl(nmDataCurr, {sDs.varDs, 'sim','hist'}, 'and') && ~regexpbl(nmDataCurr, '4bc')
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
end


%ENSURE VARIABLES TO USE IN BIAS CORRECTION MATCH THE DOWNSCALING VARIABLE
if ~isequal(varBase, sDs.varDs)
    temp = sSimBase.(varBase);
    sSimBase.(sDs.varDs) = temp;
    sSimBase = rmfield(sSimBase, varBase);
end

if ~isequal(var2Bc, sDs.varDs)
    error('bcSwitch:diffVar',['The variable chosen for bias correction is ' ...
        var2Bc ' but the downscaling variable is ' sDs.varDs '. This is not expected.'])
end


%Implement downscaling
sDataOut = stn_eQM_bc(sStnCurr, sSimBase, sSim2Bc, sDs.varDs, sDs.bctype);