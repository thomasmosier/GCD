function dataModBcOut = bc1d_month(dataHistRef, datesHistRef, dataHistMod, datesHistMod, dataMod2Bc, datesMod2Bc, type)

dataModBcOut = nan(size(dataMod2Bc), 'single');

if sum(size(dataHistRef) > 1) > 1
    error('bc1dMnth:histRefNotVec','The historic reference data is not a vector. Input must be a vector.')
end
if sum(size(dataHistMod) > 1) > 1
    error('bc1dMnth:histModNotVec','The historic model data is not a vector. Input must be a vector.')
end
if sum(size(dataMod2Bc) > 1) > 1
    error('bc1dMnth:mod2BcNotVec','The model data being bias corrected is not a vector. Input must be a vector.')
end

mnths = unique(datesMod2Bc(:,2));
for kk = 1 : numel(mnths)
    indRefHCurr = find(datesHistRef(:,2) == mnths(kk));
    indSimHCurr = find(datesHistMod(:,2) == mnths(kk));
    indSimPCurr = find( datesMod2Bc(:,2) == mnths(kk));
    
    if regexpbl(type, 'eQM')
        [dataModBcOut(indSimPCurr)] = e_qm(...
            dataHistRef(indRefHCurr), ...
            dataHistMod(indSimHCurr), ...
            dataMod2Bc( indSimPCurr), type);
    else regexpbl(type, 'q2q')
        [dataModBcOut(indSimPCurr)] = e_q2q(...
            dataHistRef(indRefHCurr), ...
            dataHistMod(indSimHCurr), ...
            dataMod2Bc( indSimPCurr));
    end
end