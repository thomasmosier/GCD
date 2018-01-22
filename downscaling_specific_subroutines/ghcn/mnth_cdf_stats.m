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

function stats = mnth_cdf_stats(gridData, stnRec, nBins, var, path)




if regexpbl(var,'pr')
    units = 'mm'; 
elseif regexpbl(var,{'tmp', 'tmn', 'tmx'})
   units = 'Celsius'; 
else
   units = 'Unknown'; 
end
[~, indSort] = sort(cell2mat(stnRec(:,4)));

gridData = gridData(indSort,:);
stnRec = stnRec(indSort,:);

nSeries = 12;
%%SET PROPERTIES:
%Set font size:
ftSz = 14;
    ftSzT = ftSz + 2;
    ftSzAx = ftSz - 2;
%Set line width:
lnWd = 2;
    lnWdA = lnWd - 0.5;
strFont = 'Arial';
%Use these colors instead?
%From color Brewer):
clrBrewer = [228,26,28; ...
    55,126,184; ...
    77,175,74; ...
    152,78,163]/255;
% szFrame = [3.5,3];
szFrame = [12,7];


mnthBias  = nan(13,1);
mnthMAE   = nan(13,1);
mnthMAPE  = nan(13,1);
mnthWMAPE = nan(13,1);
mnthCorr = nan(13,2);
nMonthPts = nan(13,1);
valGrid = cell(13,1);
valStn = cell(13,1);
% stnNse   = nan(length(gridData(:,1)),1);
% stnKge   = nan(length(gridData(:,1)),1);
    

%Set colors:
colorsUse = distinguishable_colors( nSeries );

nMnths = 12;
%Loop over all months:
for ii = 1 : nMnths
    
    mnthGridVec = [];
    mnthStnVec = [];
    for jj = 1 : length(stnRec(:,1))
        mnthGridVec = cat(1, mnthGridVec, ...
            reshape(gridData{jj,5}(:, ii + 1), [], 1));
        mnthStnVec = cat(1, mnthStnVec, ...
            reshape(stnRec{jj,5}(:, ii + 1), [], 1));
    end
    
    nBins = (min([nBins, numel(mnthGridVec), numel(mnthStnVec)]));
    [cdfGrid, valGrid{ii}] = e_cdf(mnthGridVec, 'bins', nBins);
    [cdfStn,   valStn{ii}] = e_cdf( mnthStnVec, 'bins', nBins);
    
    if all2d(isnan(valGrid{ii})) || all2d(isnan(valStn{ii}))
        continue
    end
        
    nMonthPts(ii) = sum(~isnan(mnthGridVec));
    
    if isnan(nMonthPts(ii)) || isempty(nMonthPts(ii)) || nMonthPts(ii) == 0
        mnthBias(ii)  = nan;
        mnthMAE(ii)   = nan;
        mnthMAPE(ii)  = nan;
        mnthWMAPE(ii) = nan;
        mnthCorr(ii,:) = nan(1,2);
    else
        mnthBias(ii)  = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),  'bias');
        mnthMAE(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'MAE');
        mnthMAPE(ii)  = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),  'MAPE');
        mnthWMAPE(ii) = fitness(valStn{ii}(2:end), valGrid{ii}(2:end), 'WMAPE');
        mnthCorr(ii,:) = fitness(valStn{ii}(2:end), valGrid{ii}(2:end), 'pearson');
    end
     
%     max(reshape(gridData{ii,5}(:,2:13),1,[]))
%     max(reshape(stnRec{ii,5}(:,2:13),1,[]))
    

%     stnNse(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'nse');
%     stnKge(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'kge');
    

    iterMod = mod(ii-1, nSeries)+1;

    if iterMod == 1
        strSeries = cell(2*nSeries,1);
        hTs = nan(2*nSeries,1);
        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'visible', 'off');
        hold on
    end
    
    
    indCurr = [2*(iterMod-1) + 1, 2*iterMod];
    hTs(indCurr) = plot(...
         valStn{ii}(2:end),  cdfStn(2:end), '-', ...
         valGrid{ii}(2:end), cdfGrid(2:end), ':');
     
     %Set properties of current series:
    set(hTs(indCurr),'LineWidth',lnWd, ...
        'Color',colorsUse(iterMod,:));
     
    strSeries{indCurr(1)} = [char(num2Month(ii, 3)) ' (Stn)'];
    strSeries{indCurr(2)} = [char(num2Month(ii, 3)) ' (Grid)'];
     
     
     if iterMod == nSeries || ii == nMnths % | ii == length(gridData(:,1))       
        %Create legend:
        indRemLgd = find(isnan(hTs) == 1);
        hTs(indRemLgd) = [];
        strSeries(indRemLgd) = [];
        hLgd = legend(hTs, strSeries);
        
        hYLab = ylabel('Cumulative Probability');
        hXLab = xlabel([var ' (' units ')']);
        hTtl = title('Comparison of Monthly CDFs for GHCN Records Relative to Downscaled Data');
        
        
        %Set axis and data frame properties:
        set(hLgd, ...
            'FontSize'   , ftSz, ...
            'LineWidth', lnWdA,...
            'location','SouthEastOutside');
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
            ['CDF_mnth_', var]);
        %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        savefig(hFig,[pathCdfPlot '.fig']);
        print(hFig,[pathCdfPlot '.eps'],'-depsc2');
        print(hFig,[pathCdfPlot '.png'],'-dpng','-r600');
        
        
        hold off
     end
end
hold off


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
% 
% stats(5).name = 'Stn NSE';
% stats(5).data = stnNse;
% 
% stats(6).name = 'Stn KGE';
% stats(6).data = stnKge;


