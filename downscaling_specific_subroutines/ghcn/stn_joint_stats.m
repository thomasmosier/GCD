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

function [stats, corCom] = stn_joint_stats(gridReform, stnRec, nBins, dirReport)

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


%Create figure for EACH STATION:
stnStnData = cell(2,numel(stnRec{1}(:,1)));
[stnStnData{:}] = deal([]);
stnGridData = cell(2,numel(stnRec{1}(:,1)));
[stnGridData{:}] = deal([]);

stnBias  = nan(1+length(gridReform{1}(:,1)),1);
stnRMSE  = nan(1+length(gridReform{1}(:,1)),1);
stnMAE   = nan(1+length(gridReform{1}(:,1)),1);
stnMAPE  = nan(1+length(gridReform{1}(:,1)),1);
stnWMAPE = nan(1+length(gridReform{1}(:,1)),1);
stnCorr = nan(1+length(gridReform{1}(:,1)),2);
stnCorGHCN = nan(1+length(gridReform{1}(:,1)),2);
stnCorGrid = nan(1+length(gridReform{1}(:,1)),2);


%Put data into station aggregated form
for kk = 1 : 2
    for ii = 1 : numel(stnRec{1}(:,1))
            stnStnData{kk,ii}  =     stnRec{kk}{ii,5}(:,2:13);
            stnGridData{kk,ii} = gridReform{kk}{ii,5}(:,2:13);
    end
end

for jj = 1 : numel(stnRec{1}(:,1))
    indNaNGHCN = find(isnan(stnGridData{1,jj}) | isnan(stnGridData{2,jj}));
    stnGridData{1,jj}(indNaNGHCN) = [];
    stnGridData{2,jj}(indNaNGHCN) = [];
    indNaNStn = find(isnan(stnStnData{1,jj}) | isnan(stnStnData{2,jj}));
    stnStnData{1,jj}(indNaNStn) = [];
    stnStnData{2,jj}(indNaNStn) = [];
end


for ii = 1 : numel(stnRec{1}(:,1))
    [~, val{1}] = e_cdf(stnStnData{1,ii}, 'bins', nBins);
    [~, val{2}] = e_cdf(stnStnData{2,ii}, 'bins', nBins);

    [nJointStn, ~, ~] = joint_hist(stnStnData{1,ii},stnStnData{2,ii}, val{1}, val{2});
    [nJointGrid, ~, ~] = joint_hist(stnGridData{1,ii},stnGridData{2,ii}, val{1}, val{2});

    if all2d(isnan(nJointStn)) || all2d(isnan(nJointGrid))
        continue
    end
    
    hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame],'visible','off');

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
            ['Joint_hist', '_stn_', num2str(stnRec{1}{ii,1})]);
    
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hFig,[pathHistPlot '.fig']);
    print(hFig,[pathHistPlot '.eps'],'-depsc2');
    print(hFig,[pathHistPlot '.png'],'-dpng','-r600');


    stnBias(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'bias');
    stnRMSE(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'RMSE');
    stnMAE(ii)   = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),   'MAE');
    stnMAPE(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'MAPE');
    stnWMAPE(ii) = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'WMAPE');
    stnCorr(ii,:) = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'pearson');
    
    %Write correlations:
    [rGHCN,pGHCN] = corrcoef(stnStnData{1,ii},stnStnData{2,ii});
    [rGrid,pGrid] = corrcoef(stnGridData{1,ii},stnGridData{2,ii});
    stnCorGHCN(ii,:) = [rGHCN(2), pGHCN(2)];
    stnCorGrid(ii,:) = [rGrid(2), pGrid(2)];
end

stnBias(end)   = nanmean(stnBias(1:end-1));
stnRMSE(end)   = nanmean(stnRMSE(1:end-1));
stnMAE(end)    = nanmean(stnMAE(1:end-1));
stnMAPE(end)   = nanmean(stnMAPE(1:end-1));
stnWMAPE(end)  = nanmean(stnWMAPE(1:end-1));
stnCorr(end,:) = nanmean(stnCorr(1:end-1,:),1);
stnCorGHCN(end,:) = nanmean(stnCorGHCN(1:end-1,:),1);
stnCorGrid(end,:) = nanmean(stnCorGrid(1:end-1,:),1);


%Fill output structure:
stats(1).name = 'Stn Bias';
stats(1).data = stnBias;

stats(2).name = 'Stn RMSE';
stats(2).data = stnRMSE;

stats(3).name = 'Stn MAE';
stats(3).data = stnMAE;

stats(4).name = 'Stn MAPE';
stats(4).data = 100*stnMAPE;

stats(5).name = 'Stn WMAPE';
stats(5).data = 100*stnWMAPE;

stats(6).name = 'Stn r';
stats(6).data = stnCorr(:,1);

stats(7).name = 'Stn p-val';
stats(7).data = stnCorr(:,2);

% 
% stats(5).name = 'Stn NSE';
% stats(5).data = stnNse;
% 
% stats(6).name = 'Stn KGE';
% stats(6).data = stnKge;


%Write GHCN correlation output:
corCom(1).name = 'Stn GHCN r';
corCom(1).data = stnCorGHCN(:,1);

corCom(2).name = 'Stn GHCN p-val';
corCom(2).data = stnCorGHCN(:,2);

%Write Gridd correlation output:
corCom(3).name = 'Stn Grid r';
corCom(3).data = stnCorGrid(:,1);

corCom(4).name = 'Stn Grid p-val';
corCom(4).data = stnCorGrid(:,2);

corCom(5).name = 'MAE in r';
corCom(5).data = [abs(stnCorGrid(1:end-1,1) - stnCorGHCN(1:end-1,1)); nanmean(abs(stnCorGrid(1:end-1,1) - stnCorGHCN(1:end-1,1)))];

