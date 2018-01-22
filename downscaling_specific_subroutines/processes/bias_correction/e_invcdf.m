function cdfNew = e_invcdf(dataRef, dataNew)

%Initialize output:
cdfNew = nan(size(dataNew));

[cdfRef, valRef] = e_cdf(dataRef);

for ii = 1 : numel(dataNew)
   [~, indNew] = min(abs(valRef - dataNew(ii)));
   cdfNew(ii) = cdfRef(indNew);
end

