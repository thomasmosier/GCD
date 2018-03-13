% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.

function stats = stn_stats_calc_v2(gridData, stnData)


stats = struct;

stnBias  = nan(1+length(gridData(:,1)),1);
stnMAE   = nan(1+length(gridData(:,1)),1);
stnMAPE  = nan(1+length(gridData(:,1)),1);
stnWMAPE = nan(1+length(gridData(:,1)),1);
stnCorr  = nan(1+length(gridData(:,1)),2);

for ii = 1 : length(gridData(:,1))

    stnVec = stnData{ii,5}(:,2:13);
    gridVec = gridData{ii,5}(:,2:13);
    indUse = find(~isnan(stnVec) & ~isnan(gridVec));
    
    if isempty(indUse) || numel(indUse) == 0
        continue
    else 
        stnBias(ii)   = fitness(stnVec(indUse), gridVec(indUse),  'bias');
        stnMAE(ii)    = fitness(stnVec(indUse), gridVec(indUse),   'MAE');
        stnMAPE(ii)   = fitness(stnVec(indUse), gridVec(indUse),  'MAPE');
        stnWMAPE(ii)  = fitness(stnVec(indUse), gridVec(indUse), 'WMAPE');
        stnCorr(ii,:) = fitness(stnVec(indUse), gridVec(indUse), 'pearson');
%         %BIAS
%         stnBias(ii) = nansumn(nansumn( stnTempDiff )) / nStnPts(ii);
% 
%         %MEAN ABSOLUTE ERROR (MAE)
%         stnMae(ii) = nansumn(nansumn( abs(stnTempDiff) )) / nStnPts(ii);
% 
%         %MAPE
%         stnMapeTemp = abs(stnTempDiff./stnTotTemp);
%         stnMapeNZero = sum(sum(stnMapeTemp == Inf));
%         stnMapeTemp(stnMapeTemp == Inf) = NaN;
%         stnMape(ii) = (100/(nStnPts(ii) - stnMapeNZero))*nansumn(nansumn( stnMapeTemp ));
% 
%         %RMSE
%         stnRmse(ii) = sqrt( nansumn(nansumn(stnTempDiff.^2  )) / nStnPts(ii) );
% 
%         %WEIGHTED MEAN ABSOLUTE PERCENT ERROR
%          stnWmape(ii) = 100 * nansumn(nansumn( abs( stnTempDiff ) )) ...
%              / nansumn(nansumn( abs(stnTotTemp) ));    
% 
%         %NSE
%         stnNseAvg = stnTotTemp / nStnPts(ii);
%         stnNse(ii) = 1 ...
%             - nansumn(nansumn( stnTempDiff.^2 )) ...
%             / nansumn(nansumn( (stnTotTemp - stnNseAvg).^2 ));
    end
end

stnBias(end)   = nanmean(stnBias(1:end-1));
stnMAE(end)    = nanmean(stnMAE(1:end-1));
stnMAPE(end)   = nanmean(stnMAPE(1:end-1));
stnWMAPE(end)  = nanmean(stnWMAPE(1:end-1));
stnCorr(end,:) = nanmean(stnCorr(1:end-1,:),1);

%Fill output structure:
stats(1).name = 'Stn Bias';
stats(1).data = stnBias;

stats(2).name = 'Stn MAE';
stats(2).data = stnMAE;

stats(3).name = 'Stn MAPE';
stats(3).data = 100*stnMAPE;

stats(4).name = 'Stn WMAPE';
stats(4).data = 100*stnWMAPE;

stats(5).name = 'Stn r';
stats(5).data = stnCorr(:,1);

stats(6).name = 'Stn p-val';
stats(6).data = stnCorr(:,2);

end
