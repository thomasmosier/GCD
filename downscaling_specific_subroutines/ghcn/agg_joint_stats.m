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

function [stats, corCom] = agg_joint_stats(gridReform, stnRec, nBins, dirReport)

%Create figure for all region data AGGREGATED:
 %%SET PROPERTIES:
%Set font size:
ftSz = 14;
    ftSzT = ftSz + 2;
    ftSzAx = ftSz - 2;
%Set line width:
lnWd = 2;
    lnWdA = lnWd - 0.5;
strFont = 'Arial';
szFrame = [7,4];
tickPts = linspace(1,nBins,10);

%Create aggregated jointed probability distribution:
aggStnData = cell(2,1);
[aggStnData{:}] = deal([]);
aggGridData = cell(2,1);
[aggGridData{:}] = deal([]);

%Put data into aggregated form
    for ii = 1 : numel(stnRec{1}(:,1))
        aggStnData{1}  = cat(1, aggStnData{1}, reshape(    stnRec{1}{ii,5}(:,2:13),[],1));
        aggStnData{2}  = cat(1, aggStnData{2}, reshape(    stnRec{2}{ii,5}(:,2:13),[],1));
        aggGridData{1} = cat(1, aggGridData{1}, reshape(gridReform{1}{ii,5}(:,2:13),[],1));
        aggGridData{2} = cat(1, aggGridData{2}, reshape(gridReform{2}{ii,5}(:,2:13),[],1));
    end

indNaNGrid = find(isnan(aggGridData{1}) | isnan(aggGridData{2}));
aggGridData{1}(indNaNGrid) = [];
aggGridData{2}(indNaNGrid) = [];
indNaNGHCN = find(isnan(aggStnData{1}) | isnan(aggStnData{2}));
aggStnData{1}(indNaNGHCN) = [];
aggStnData{2}(indNaNGHCN) = [];

[~, val{1}] = e_cdf(aggStnData{1}, 'bins', nBins);
[~, val{2}] = e_cdf(aggStnData{2}, 'bins', nBins);

[ nJointStn, ~, ~] = joint_hist( aggStnData{1}, aggStnData{2}, val{1}, val{2});
[nJointGrid, ~, ~] = joint_hist(aggGridData{1},aggGridData{2}, val{1}, val{2});

hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Visible', 'off');

hSub1 = subplot(1,2,1);
pcolor((1:nBins), (1:nBins), nJointStn);shading flat
hXLab{2} = xlabel('Temperature (Deg C)'); hYLab{1} = ylabel('Precipitation (mm)');hTtl{1} = title('GHCN Station data');
ax = gca;
ax.XTick = tickPts;
ax.YTick = tickPts;
ax.XTickLabel = round2(val{1}(round(tickPts)),0);
ax.YTickLabel = round2(val{2}(round(tickPts)),0);

hSub2 = subplot(1,2,2);
pcolor((1:nBins), (1:nBins), nJointGrid);shading flat
hXLab{2} = xlabel('Temperature (Deg C)'); hYLab{2} = ylabel('Precipitation (mm)'); hTtl{2} = title('Gridded Data');
ax = gca;
ax.XTick = tickPts;
ax.YTick = tickPts;
ax.XTickLabel = round2(val{1}(round(tickPts)),0);
ax.YTickLabel = round2(val{2}(round(tickPts)),0);

s1Pos = get(hSub1,'position'); % [left bottom width height]
s2Pos = get(hSub2,'position'); % [left bottom width height]

hCbAx = axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
hCb = colorbar ('location','southoutside','FontSize',ftSz);
maxP =  max(max2d(nJointStn),max2d(nJointGrid));
caxis(hSub1, [0, maxP]);
caxis(hSub2, [0, maxP]);
caxis(hCbAx, [0, maxP]);
cbT = round2(linspace(0, maxP, 11),1);
set(hCb,'TickLabels',cbT);

s1Pos(4) = s1Pos(4) - s1Pos(3);
s1Pos(2) = s1Pos(2) +( s1Pos(4) - s1Pos(3)) + 0.1;

s2Pos(2) = s1Pos(2);
s2Pos(4) = s1Pos(4);
s2Pos(1) = s2Pos(1) + 0.05;

set(hSub1,'position',s1Pos);
set(hSub2,'position',s2Pos);

ylabel(hCb,'Joint Probability Density (%)');
cpos = hCb.Position;
cpos(4) = 0.5*cpos(4);
cpos(2) = cpos(2) + 0.05;
hCb.Position = cpos;
%     hCb = colorbar(gca,'location','southoutside');

set([hXLab{:}, hYLab{:}]  , ...
    'FontSize'   , ftSz, ...
    'Color', 'black', ...
    'FontName'   , strFont);
set([hTtl{:}], ...
    'FontWeight' , 'bold', ...
    'FontSize'   , ftSzT);
set([hSub1, hSub2], ...
    'Box'         , 'on', ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'YMinorTick'  , 'off'      , ...
    'fontSize'    , ftSzAx, ...
    'LineWidth'   , lnWdA, ...
    'FontName'   , strFont);
set(hCb, ...
    'fontSize'    , ftSzAx, ...
    'LineWidth'   , lnWdA, ...
    'FontName'   , strFont);


%Export:
pathHistPlot = fullfile(dirReport, ...
    'Joint_hist_agg');
%set(gcf,'units','normalized','outerposition',[0 0 1 1]);
savefig(hFig,[pathHistPlot '.fig']);
print(hFig,[pathHistPlot '.eps'],'-depsc2');
print(hFig,[pathHistPlot '.png'],'-dpng','-r600');


aggBias  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'bias');
aggRMSE  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'RMSE');
aggMAE   = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),   'MAE');
aggMAPE  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'MAPE');
aggWMAPE = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'WMAPE');
aggCorr  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'pearson');


%Fill output structure:
stats(1).name = 'Agg Bias';
stats(1).data = aggBias;

stats(2).name = 'Agg RMSE';
stats(2).data = aggRMSE;

stats(3).name = 'Agg MAE';
stats(3).data = aggMAE;

stats(4).name = 'Agg MAPE';
stats(4).data = 100*aggMAPE;

stats(5).name = 'Agg WMAPE';
stats(5).data = 100*aggWMAPE;

stats(6).name = 'Agg r';
stats(6).data = aggCorr(1);

stats(7).name = 'Agg p-val';
stats(7).data = aggCorr(2);

[rGHCN,pGHCN] = corrcoef(aggStnData{1},aggStnData{2});
corCom(1).name = 'Agg GHCN r';
corCom(1).data = rGHCN(2);

corCom(2).name = 'Agg GHCN p-val';
corCom(2).data = pGHCN(2);

[rGrid,pGrid] = corrcoef(aggGridData{1},aggGridData{2});
corCom(3).name = 'Agg Grid r';
corCom(3).data = rGrid(2);

corCom(4).name = 'Agg Grid p-val';
corCom(4).data = pGrid(2);

corCom(5).name = 'MAE in r';
corCom(5).data = nanmean(abs(rGrid(2)-rGHCN(2)));

