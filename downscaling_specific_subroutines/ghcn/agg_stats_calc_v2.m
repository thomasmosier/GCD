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

aggStats = struct;
gridVec = [];
stnVec = [];

for ii = 1 : length(gridData(:,1))
    gridVec = [gridVec; reshape(gridData{ii,5}(:,2:13),[], 1)];
    stnVec  = [ stnVec; reshape( stnData{ii,5}(:,2:13),[], 1)];
end
clear ii

indUse = find(~isnan(gridVec) & ~isnan(stnVec));

if numel(indUse) == 0
    totBias  = nan;
    totMae   = nan;
    totRmse  = nan;
    totNse   = nan;
    totMape  = nan;
    totWmape = nan;
    totCorr  = nan(1,2);
else
    gridVec = gridVec(indUse);
    stnVec  = stnVec(indUse);

    totBias  = fitness(stnVec, gridVec,  'bias');
    totMae   = fitness(stnVec, gridVec,   'MAE');
    totRmse  = fitness(stnVec, gridVec,  'rmse');
    totNse   = fitness(stnVec, gridVec,  'nse');
    totMape  = fitness(stnVec, gridVec,  'MAPE');
    totWmape = fitness(stnVec, gridVec, 'WMAPE');
    totCorr  = fitness(stnVec, gridVec, 'pearson');
end

%Fill output structure:
aggStats(1).name = 'Agg Bias';
aggStats(1).data = totBias;

aggStats(2).name = 'Agg MAE';
aggStats(2).data = totMae;

aggStats(3).name = 'Agg RMSE';
aggStats(3).data = totRmse;

aggStats(4).name = 'Agg NSE';
aggStats(4).data = 1 - totNse;

aggStats(5).name = 'Agg MAPE';
aggStats(5).data = 100*totMape;

aggStats(6).name = 'Agg WMAPE';
aggStats(6).data = 100*totWmape;

stats(6).name = 'Agg r';
stats(6).data = totCorr(1);

stats(6).name = 'Agg p-val';
stats(6).data = totCorr(2);



