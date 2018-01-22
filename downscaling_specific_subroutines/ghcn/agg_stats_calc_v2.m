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

function aggStats = agg_stats_calc_v2(gridData, stnData)


warning('agg_stats_cal_v2:outdated','This function must be upgraded for it to work properly.')
    nTotPts = 0;
totBias = 0;
totMape = 0;
    totMapeNZero = 0;
totRmse = 0;
totWmape = 0;
    totStnSum = 0;
totNse = 0;
    totNseDenTemp = [];
totMae = 0;

gridVec = [];
stnVec = [];
    
for ii = 1 : length(gridData(:,1))
    gridVec = cat(1, gridVec, reshape(gridData{ii,5}(:,2:13),1,[]));
    stnVec = cat(1, stnVec, reshape(stnData{ii,5}(:,2:13),1,[]));
end

totTempDiff = gridVec - stnVec;

nTotPts = sum(sum(~isnan(stnVec)));

%MEAN ERROR (BIAS)
totBias = nansumn(nansumn( totTempDiff , 1));

%MEAN ABSOLUTE ERROR (MAE)
totMae = nansumn(nansumn( abs(totTempDiff) ));

%Mean ABSOLUTE PERCENT ERROR (MAPE)
totMapeTemp = abs(totTempDiff ./ stnVec );
totMapeIndZero = find(totMapeTemp == Inf);  %These points occur where station measurements are exactly 0.
totMapeNZero = totMapeNZero + length(totMapeIndZero);
totMapeTemp( totMapeTemp == Inf ) = NaN;
totMape = totMape + nansumn(nansumn( totMapeTemp ));

%ROOT MEAN SQUARE ERROR (RMSE)
totRmse = totRmse + nansumn(nansumn(totTempDiff.^2));

%WEIGHTED MEAN ABSOLUTE PERCENT ERROR (WMAPE)
totWmape = totWmape + nansumn(nansumn( abs( totTempDiff ) ));

totStnSumTemp = stnData{ii,5}(:,2:13);
totStnSumTemp = totStnSumTemp(~isnan(totTempDiff));
totStnSum = nansumn(nansumn( abs(totStnSumTemp) ));

%NASH-SUTCLIFFE EFFICIENCY (NSE)
totNse = nansumn(nansumn( totTempDiff.^2 ));
totNseDenTemp = cat(1, totNseDenTemp, reshape(totStnSumTemp, [], 1) );

if isnan(nTotPts) || isempty(nTotPts) || nTotPts == 0
    totBias  = nan;
    totMae   = nan;
    totRmse  = nan;
    totMape  = nan;
    totWmape = nan;
    totNse   = nan;
else
    totBias = totBias / nTotPts;
    totBias(totBias == 0) = NaN;
    
    totMae = totMae / nTotPts;
    totMae(totMae == 0) = NaN;
    
    totMape = (100/(nTotPts - totMapeNZero)) * totMape;
    totMape(totMape == 0) = NaN;
    
    totRmse = sqrt( totRmse / nTotPts );
    totRmse(totRmse == 0) = NaN;

    totWmape = 100 * totWmape / totStnSum;
    totWmape(totWmape == 0) = NaN;

    totStnAvg = totStnSum / nTotPts;
    totNseDen = nansumn( (totNseDenTemp - totStnAvg).^2 );
    totNse = 1 - totNse / totNseDen;
    totNse(totNse == 0) = NaN;
end

%Fill output structure:
aggStats(1).name = 'Agg Bias';
aggStats(1).data = totBias;

aggStats(2).name = 'Agg MAE';
aggStats(2).data = totMae;

aggStats(3).name = 'Agg RMSE';
aggStats(3).data = totRmse;

aggStats(4).name = 'Agg NSE';
aggStats(4).data = totNse;

aggStats(5).name = 'Agg MAPE';
aggStats(5).data = 100*totMape;

aggStats(6).name = 'Agg WMAPE';
aggStats(6).data = 100*totWmape;




