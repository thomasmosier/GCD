function [dataModCorr] = e_q2q(dataHistRef, dataHistMod, dataMod2Bc)


%Initialize output:
dataModCorr = nan(size(dataMod2Bc));

if all(isnan(dataHistRef)) || all(isnan(dataHistMod))
    dataModCorr = nan(size(dataMod2Bc));
    return
end


[cdfHistRef, valHistRef] = e_cdf(dataHistRef, 'bins', 100);
[cdfHistMod, valHistMod] = e_cdf(dataHistMod, 'bins', 100);

%COMPARE BINNED AND UNBINNED CDFs
% [cdfHistRefNB, valHistRefNB] = e_cdf(dataHistRef);
% [cdfHistModNB, valHistModNB] = e_cdf(dataHistMod);
% figure; hPlot = plot(valHistRef, cdfHistRef, '-.', valHistRefNB, cdfHistRefNB, '-', valHistMod, cdfHistMod, '-.', valHistModNB, cdfHistModNB, '-', 'linewidth', 1.5); legend(hPlot, {'Historic Ref (binned)', 'Historic Ref', 'Historic Sim (binned)', 'Historic Sim'}, 'location', 'southeast');
% xlim([min([valHistRef(:); valHistMod(:); valHistRefNB(:); valHistModNB(:)]), max([valHistRef(:); valHistMod(:); valHistRefNB(:); valHistModNB(:)])]);


%COMPARE ORIGINAL CDF TO KERNPDF CDF
% valHistRef = sort(100*rand(101,1));
% valHistMod = sort(10*rand(101,1));

% valHistRefKern = kernpdf_modify(valHistRef, cdfHistMod);
% valHistModKern = kernpdf_modify(valHistMod, cdfHistMod);
% figure; hPlot = plot(valHistMod, cdfHistMod, valHistModKern, cdfHistMod, 'linewidth', 1.5); legend(hPlot, {'Original CDF', 'Kernel CDF'},'location', 'southeast');
% figure; hPlot = plot(valHistRef, cdfHistRef, valHistRefKern, cdfHistRef, 'linewidth', 1.5); legend(hPlot, {'Original CDF', 'Kernel CDF'},'location', 'southeast');


minHistMod = min(valHistMod);
maxHistMod = max(valHistMod);

for kk = 1 : numel(dataMod2Bc)
    if isnan(dataMod2Bc(kk))
        dataModCorr(kk) = nan;
    else
        if dataMod2Bc(kk) >= maxHistMod %Outside of range (Use highest quantile)
            bias = valHistMod(end) - valHistRef(end);
        elseif dataMod2Bc(kk) <= minHistMod %Outside of range (Use lowest quantile)
            bias = valHistMod(1) - valHistRef(1);
        else %Inside of range
            indModGT = find(dataMod2Bc(kk) < valHistMod, 1, 'first');
            
            if indModGT == 1 %This case shouldn't occur because it would be same as outside of range (use lowest quantile)
                bias = valHistMod(1) - valHistRef(1);
            else %Interpolate between neighboring quantiles
                indModLT = indModGT - 1;
            
                %Estimate cdf value of current data based on linear weighting of historic values:
                wgtModLT = (valHistMod(indModGT) - dataMod2Bc(kk)) / (valHistMod(indModGT) - valHistMod(indModLT));
                wgtModGT = (dataMod2Bc(kk) - valHistMod(indModLT)) / (valHistMod(indModGT) - valHistMod(indModLT));
                
                if round2(wgtModLT+wgtModGT, 4) ~= 1
                   error('q2q:wgtsWrong',['The weights sum to ' num2str(wgtModLT+wgtModGT) ' instead of 1.']); 
                end
                
                valModCurr = wgtModLT*valHistMod(indModLT) + wgtModGT*valHistMod(indModGT);
                cdfMod2BcCurr = wgtModLT*cdfHistMod(indModLT) + wgtModGT*cdfHistMod(indModGT);
                
                %Find corresponding CDF in reference data:
                indRefGT = find(cdfMod2BcCurr < cdfHistRef, 1, 'first');
                
                if isempty(indRefGT) 
                    valRefCurr = valHistRef(end);
                elseif indRefGT == 1
                    valRefCurr = valHistRef(1);
                else %Linear weight reference data:
                    indRefLT = indRefGT - 1;
                    wgtRefLT = (cdfHistRef(indRefGT) - cdfMod2BcCurr) / (cdfHistRef(indRefGT) - cdfHistRef(indRefLT));
                    wgtRefGT = (cdfMod2BcCurr - cdfHistRef(indRefLT)) / (cdfHistRef(indRefGT) - cdfHistRef(indRefLT));
                    
                    if round2(wgtRefLT+wgtRefGT, 4) ~= 1
                        error('q2q:wgtsWrong',['The weights sum to ' num2str(wgtRefLT+wgtRefGT) ' instead of 1.']); 
                    end
                
                    valRefCurr = wgtRefLT*valHistRef(indRefLT) + wgtRefGT*valHistRef(indRefGT);
                end

                %Bias correction is a weighted mean of the two points
                bias = valModCurr - valRefCurr;
            end
        end
        
        dataModCorr(kk) = dataMod2Bc(kk) - bias; 
    end
end


%COMPARE INPUT AND OUTPUT CDFS
% [cdfMod2Bc, valMod2Bc] = e_cdf(dataMod2Bc, 'bins', 100);
% [cdfModCorr, valModCorr] = e_cdf(dataModCorr, 'bins', 100);
% figure; hPlot = plot(valHistRef, cdfHistRef, '-', valHistMod, cdfHistMod, '-', valMod2Bc, cdfMod2Bc, '.', valModCorr, cdfModCorr, '*', 'linewidth', 1.5); legend(hPlot, {'Historic Ref', 'Historic Sim', 'Sim 2 BC', 'Sim Correct'}, 'location', 'southeast');
% xlim([min([valHistRef(:); valHistMod(:); valMod2Bc(:); valModCorr(:)]), max([valHistRef(:); valHistMod(:); valMod2Bc(:); valModCorr(:)])]);


