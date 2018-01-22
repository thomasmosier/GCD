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

function stats = agg_cdf_stats(gridData, stnData, nBins, var, path)


if regexpbl(var,'pr')
    units = 'mm'; 
elseif regexpbl(var,{'tmp', 'tmn', 'tmx'})
   units = 'Celsius'; 
else
   units = 'Unknown'; 
end

%%SET PROPERTIES:
%Set font size:
ftSz = 14;
    ftSzT = ftSz + 2;
    ftSzAx = ftSz - 2;
%Set line width:
lnWd = 2;
    lnWdA = lnWd - 0.5;
strFont = 'Arial';

% szFrame = [3.5,3];
szFrame = [12,7];
    
gridVec = [];
stnVec = [];

%Loop over all stations:
for ii = 1 : length(gridData(:,1))
    gridVec = cat(1, gridVec, reshape(gridData{ii,5}(:,2:13),[],1));
    stnVec = cat(1, stnVec, reshape(stnData{ii,5}(:,2:13),[],1));
end

[cdfGrid, valGrid] = e_cdf(gridVec, 'bins', nBins);
[cdfStn, valStn]  = e_cdf(stnVec, 'bins', nBins);

aggBias  = fitness(valStn(2:end), valGrid(2:end),  'bias');
aggMAE   = fitness(valStn(2:end), valGrid(2:end),   'MAE');
aggMAPE  = fitness(valStn(2:end), valGrid(2:end),  'MAPE');
aggWMAPE = fitness(valStn(2:end), valGrid(2:end), 'WMAPE');
aggCorr = fitness(valStn(2:end), valGrid(2:end), 'pearson');



%     stnNse(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'nse');
%     stnKge(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'kge');

hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'visible', 'off');

hTs = plot(...
     valStn(2:end),  cdfStn(2:end), '-', ...
     valGrid(2:end), cdfGrid(2:end), ':');

%Set properties of current series:
set(hTs,'LineWidth',lnWd);

%Create legend:
hLgd = legend(hTs, {'Station', 'Grid'});

hYLab = ylabel('Cumulative Probability');
hXLab = xlabel([var ' (' units ')']);
hTtl = title('Comparison of Aggregatecd CDFs for GHCN Records Relative to Downscaled Data');


%Set axis and data frame properties:
set(hLgd, ...
    'FontSize'   , ftSz, ...
    'LineWidth', lnWdA,...
    'location','SouthEast');
%         set(hMnthELgd,'Layer','top');
set([hXLab, hYLab]  , ...
    'FontSize'   , ftSz, ...
    'Color', 'black', ...
    'FontName'   , strFont);
set(hTtl, ...
    'FontWeight' , 'bold', ...
    'FontSize'   , ftSzT);
set(gca, ...
    'Box'         , 'off', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'YMinorTick'  , 'on'      , ...
    'fontSize'    , ftSzAx, ...
    'LineWidth'   , lnWdA, ...
    'FontName'   , strFont);

if ~exist(path,'dir')
    mkdir(path)
end

%Export:
pathCdfPlot = fullfile(path, ...
    ['CDF_agg_', var]);
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savefig(hFig,[pathCdfPlot '.fig']);
print(hFig,[pathCdfPlot '.eps'],'-depsc2');
print(hFig,[pathCdfPlot '.png'],'-dpng','-r600');




%Fill output structure:
stats(1).name = 'Agg Bias';
stats(1).data = aggBias;

stats(2).name = 'Agg MAE';
stats(2).data = aggMAE;

stats(3).name = 'Agg MAPE';
stats(3).data = 100*aggMAPE;

stats(4).name = 'Agg WMAPE';
stats(4).data = 100*aggWMAPE;

stats(5).name = 'Agg r';
stats(5).data = aggCorr(1);

stats(6).name = 'Agg p-val';
stats(6).data = aggCorr(2);

 
% 
% stats(5).name = 'Stn NSE';
% stats(5).data = stnNse;
% 
% stats(6).name = 'Stn KGE';
% stats(6).data = stnKge;