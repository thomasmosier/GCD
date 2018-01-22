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

function [stats, corCom] = mnth_joint_stats(gridReform, stnRec, nBins, dirReport)

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

%Initialize:
mnthBias  = nan(13,1);
mnthRMSE  = nan(13,1);
mnthMAE   = nan(13,1);
mnthMAPE  = nan(13,1);
mnthWMAPE = nan(13,1);
mnthCorr = nan(13,2);
mnthCorGHCN = nan(13,2);
mnthCorGrid = nan(13,2);

% %Create aggregated jointed probability distribution:
% aggStnData = cell(2,1);
% [aggStnData{:}] = deal([]);
% aggGridData = cell(2,1);
% [aggGridData{:}] = deal([]);

%Create figure for all region BY MONTH:
mnthStnData = cell(2,12);
[mnthStnData{:}] = deal([]);
mnthGridData = cell(2,12);
[mnthGridData{:}] = deal([]);

%Put data into monthly aggregated form
for kk = 1 : 2
    for ii = 1 : numel(stnRec{1}(:,1))
        for jj = 1 : 12
            mnthStnData{kk,jj}  = cat(1,  mnthStnData{kk,jj}, reshape(    stnRec{kk}{ii,5}(:,jj+1),[],1));
            mnthGridData{kk,jj} = cat(1, mnthGridData{kk,jj}, reshape(gridReform{kk}{ii,5}(:,jj+1),[],1));
        end
    end
end

for jj = 1 : 12
    indNaNGrid = find(isnan(mnthGridData{1,jj}) | isnan(mnthGridData{2,jj}));
    mnthGridData{1,jj}(indNaNGrid) = [];
    mnthGridData{2,jj}(indNaNGrid) = [];
    indNaNGHCN = find(isnan(mnthStnData{1,jj}) | isnan(mnthStnData{2,jj}));
    mnthStnData{1,jj}(indNaNGHCN) = [];
    mnthStnData{2,jj}(indNaNGHCN) = [];
end

for ii = 1 : 12
    [~, val{1}] = e_cdf(mnthStnData{1,ii}, 'bins', nBins);
    [~, val{2}] = e_cdf(mnthStnData{2,ii}, 'bins', nBins);
    
    [nJointStn, ~, ~]  = joint_hist( mnthStnData{1,ii},mnthStnData{2,ii},  val{1}, val{2});
    [nJointGrid, ~, ~] = joint_hist(mnthGridData{1,ii},mnthGridData{2,ii}, val{1}, val{2});

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
        ['Joint_hist_mnth_' num2str(ii)]);
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hFig,[pathHistPlot '.fig']);
    print(hFig,[pathHistPlot '.eps'],'-depsc2');
    print(hFig,[pathHistPlot '.png'],'-dpng','-r600');

    mnthBias(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'bias');
    mnthRMSE(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'RMSE');
    mnthMAE(ii)   = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),   'MAE');
    mnthMAPE(ii)  = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1),  'MAPE');
    mnthWMAPE(ii) = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'WMAPE');
    mnthCorr(ii,:) = fitness(reshape(nJointStn,[],1), reshape(nJointGrid,[],1), 'pearson');
    
    %Write correlations:
    [rGHCN,pGHCN] = corrcoef(mnthStnData{1,ii},mnthStnData{2,ii});
    [rGrid,pGrid] = corrcoef(mnthGridData{1,ii},mnthGridData{2,ii});
    mnthCorGHCN(ii,:) = [rGHCN(2), pGHCN(2)];
    mnthCorGrid(ii,:) = [rGrid(2), pGrid(2)];
end


mnthBias(13)   = nanmean(mnthBias(1:end-1));
mnthRMSE(13)   = nanmean(mnthRMSE(1:end-1));
mnthMAE(13)    = nanmean(mnthMAE(1:end-1));
mnthMAPE(13)   = nanmean(mnthMAPE(1:end-1));
mnthWMAPE(13)  = nanmean(mnthWMAPE(1:end-1));
mnthCorr(13,:) = nanmean(mnthCorr(1:end-1,:),1);
mnthCorGHCN(13,:) = nanmean(mnthCorGHCN(1:end-1,:),1); 
mnthCorGrid(13,:) = nanmean(mnthCorGrid(1:end-1,:),1); 


%Fill output structure:
stats(1).name = 'Mnth Bias';
stats(1).data = mnthBias;

stats(2).name = 'Agg RMSE';
stats(2).data = mnthRMSE;

stats(3).name = 'Mnth MAE';
stats(3).data = mnthMAE;

stats(4).name = 'Mnth MAPE';
stats(4).data = 100*mnthMAPE;

stats(5).name = 'Mnth WMAPE';
stats(5).data = 100*mnthWMAPE;

stats(6).name = 'Mnth r';
stats(6).data = mnthCorr(:,1);

stats(7).name = 'Mnth p-val';
stats(7).data = mnthCorr(:,2);

%Write GHCN correlation output:
corCom(1).name = 'Mnth GHCN r';
corCom(1).data = mnthCorGHCN(:,1);

corCom(2).name = 'Mnth GHCN p-val';
corCom(2).data = mnthCorGHCN(:,2);

%Write Gridd correlation output:
corCom(3).name = 'Mnth Grid r';
corCom(3).data = mnthCorGrid(:,1);

corCom(4).name = 'Mnth Grid p-val';
corCom(4).data = mnthCorGrid(:,2);

corCom(5).name = 'MAE in r';
corCom(5).data = [abs(mnthCorGrid(1:end-1,1) - mnthCorGHCN(1:end-1,1)); nanmean(abs(mnthCorGrid(1:end-1,1) - mnthCorGHCN(1:end-1,1)))];