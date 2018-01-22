function [dataModCorr] = e_qm(dataHistRef, dataHistMod, dataMod2Bc, type)

%Initialize output:
dataModCorr = nan(size(dataMod2Bc));

if all(isnan(dataHistRef)) || all(isnan(dataHistMod))
    dataModCorr = nan(size(dataMod2Bc));
    return
end

[cdfHistRef, valHistRef] = e_cdf(dataHistRef);
[cdfHistMod, valHistMod] = e_cdf(dataHistMod);

cdfBc = e_invcdf(dataHistMod, dataMod2Bc);

limit = 3; %Limit for bias correction factor

for kk = 1 : numel(dataMod2Bc)
   [~, indHistRef] = min(abs(cdfHistRef - cdfBc(kk)));
   [~, indHistMod] = min(abs(cdfHistMod - cdfBc(kk)));
   
    if regexpbl(type, 'add')
        dataModCorr(kk) = dataMod2Bc(kk) + valHistRef(indHistRef) - valHistMod(indHistMod);
    elseif regexpbl(type, 'mult')
        bcFactor = valHistRef(indHistRef)/valHistMod(indHistMod);
        
        bcFactor = max([min([bcFactor, limit]), -limit]);
        dataModCorr(kk) = dataMod2Bc(kk)*bcFactor;
        
        if isnan(dataModCorr(kk)) && dataMod2Bc(kk) == 0 || valHistMod(indHistMod) == 0
           dataModCorr(kk) = 0;
        elseif isinf(dataModCorr(kk))
            warning('e_qm:bcInf','Infinity value created in bias correction.')
        end
    else
        error('eQm:unknownType',['Type is ' type ', but two options are mult and add.']);
    end
end


