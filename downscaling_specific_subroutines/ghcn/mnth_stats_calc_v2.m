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

function stats = mnth_stats_calc_v2(gridData, stnData)



mnthBias  = nan(13,1);
mnthMAE   = nan(13,1);
mnthMAPE  = nan(13,1);
mnthWMAPE = nan(13,1);
mnthCorr = nan(13,2);
nMonthPts = nan(13,1);

for ii = 1 : 12
    
    gridVec = nan(length(stnData(:,1)),1);
    mnthTempStnVec = nan(length(stnData(:,1)),1);
    for jj = 1 : length(stnData(:,1))
        gridVec = gridData{jj,5}(:, ii + 1);
        stnVec = stnData{jj,5}(:, ii + 1);
    end
    
    nMonthPts(ii) = sum(~isnan(mnthTempDiffVec));
    
    if isnan(nMonthPts(ii)) || isempty(nMonthPts(ii)) || nMonthPts(ii) == 0
        mnthBias(ii)  = nan;
        mnthMAE(ii)   = nan;
        mnthMAPE(ii)  = nan;
        mnthWMAPE(ii) = nan;
        mnthCorr(ii,:) = nan(1,2);
    else
        mnthBias(ii)  = fitness(stnVec{ii}(2:end), gridVec{ii}(2:end),  'bias');
        mnthMAE(ii)   = fitness(stnVec{ii}(2:end), gridVec{ii}(2:end),   'MAE');
        mnthMAPE(ii)  = fitness(stnVec{ii}(2:end), gridVec{ii}(2:end),  'MAPE');
        mnthWMAPE(ii) = fitness(stnVec{ii}(2:end), gridVec{ii}(2:end), 'WMAPE');
        mnthCorr(ii,:) = fitness(stnVec{ii}(2:end), gridVec{ii}(2:end), 'pearson');
%         %BIAS
%         mnthBias(ii) = nansumn(nansumn( mnthTempDiffVec )) / nMonthPts(ii);
% 
%         %Mean Absolute Error (MAE)
%         mnthMae(ii) = nansumn(nansumn( abs(mnthTempDiffVec) )) / nMonthPts(ii);
% 
%         %MAPE
%         mnthMapeTemp = abs(mnthTempDiffVec ./ mnthTempStnVec );
%         mnthMapeNZero = sum(mnthMapeTemp == Inf);
%         mnthMapeTemp(mnthMapeTemp == Inf) = NaN;
%         mnthMape(ii) = (100/(nMonthPts(ii) - mnthMapeNZero))*nansumn(nansumn( mnthMapeTemp ));
% 
%         %RMSE
%         mnthRmse(ii) = sqrt( nansumn(nansumn( mnthTempDiffVec.^2 )) / nMonthPts(ii) );
% 
%         %WMAPE
%          mnthWmape(ii) = 100 * nansumn(nansumn( abs(mnthTempDiffVec) )) ...
%              / nansumn(nansumn( abs(mnthTempStnVec) ));    
% 
%         %NSE
%         mnthNseAvg = mnthTempStnVec / nMonthPts(ii);
%         mnthNse(ii) = 1 - nansumn(nansumn( mnthTempDiffVec.^2 )) / nansumn(nansumn( (mnthTempStnVec - mnthNseAvg).^2 ));
    end
end



mnthBias(13)   = nanmean(mnthBias(1:end-1));
mnthMAE(13)    = nanmean(mnthMAE(1:end-1));
mnthMAPE(13)   = nanmean(mnthMAPE(1:end-1));
mnthWMAPE(13)  = nanmean(mnthWMAPE(1:end-1));
mnthCorr(13,:) = nanmean(mnthCorr(1:end-1,:),1);


%Fill output structure:
stats(1).name = 'Mnth Bias';
stats(1).data = mnthBias;

stats(2).name = 'Mnth MAE';
stats(2).data = mnthMAE;

stats(3).name = 'Mnth MAPE';
stats(3).data = 100*mnthMAPE;

stats(4).name = 'Mnth WMAPE';
stats(4).data = 100*mnthWMAPE;

stats(5).name = 'Stn r';
stats(5).data = mnthCorr(:,1);

stats(6).name = 'Stn p-val';
stats(6).data = mnthCorr(:,2);
