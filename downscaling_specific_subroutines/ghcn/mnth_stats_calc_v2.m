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

stats = struct;

mnthBias  = nan(13,1);
mnthMAE   = nan(13,1);
mnthMAPE  = nan(13,1);
mnthWMAPE = nan(13,1);
mnthRmse  = nan(13,1);
mnthCorr  = nan(13,2);

for ii = 1 : 12
    
    gridVec = [];
    stnVec  = [];
    for jj = 1 : length(stnData(:,1))
        indUse = find(~isnan(gridData{jj,5}(:, ii + 1)) & ~isnan(stnData{jj,5}(:, ii + 1)));
        gridVec = [gridVec(:); gridData{jj,5}(indUse, ii + 1)];
        stnVec  = [stnVec(:);  stnData{ jj,5}(indUse, ii + 1)];
    end

    if numel(stnVec) == 0
    	mnthBias(ii)  = nan;
        mnthMAE(ii)   = nan;
        mnthMAPE(ii)  = nan;
        mnthWMAPE(ii) = nan;
        mnthRmse(ii) = nan;
        mnthCorr(ii,:) = nan(1,2);
    else
        mnthBias(ii)   = fitness(stnVec, gridVec,  'bias');
        mnthMAE(ii)    = fitness(stnVec, gridVec,   'MAE');
        mnthMAPE(ii)   = fitness(stnVec, gridVec,  'MAPE');
        mnthWMAPE(ii)  = fitness(stnVec, gridVec, 'WMAPE');
        mnthRmse(ii)   = fitness(stnVec, gridVec, 'rmse');
        mnthCorr(ii,:) = fitness(stnVec, gridVec, 'pearson');
    end
end

%Calculate annual average values
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
