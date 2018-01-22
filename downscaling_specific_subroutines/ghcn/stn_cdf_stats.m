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

function stats = stn_cdf_stats(gridData, stnRec, nBins, var, path)



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

nSeries = 10;
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


stnBias  = nan(1+length(gridData(:,1)),1);
stnMAE   = nan(1+length(gridData(:,1)),1);
stnMAPE  = nan(1+length(gridData(:,1)),1);
stnWMAPE = nan(1+length(gridData(:,1)),1);
stnCorr = nan(1+length(gridData(:,1)),2);
valGrid = cell(1+length(gridData(:,1)),1);
valStn = cell(1+length(gridData(:,1)),1);
cntr = 0;
% stnNse   = nan(length(gridData(:,1)),1);
% stnKge   = nan(length(gridData(:,1)),1);
    

%Set colors:
colorsUse = distinguishable_colors( nSeries );


%Loop over all stations:
for ii = 1 : length(gridData(:,1))
%     max(reshape(gridData{ii,5}(:,2:13),1,[]))
%     max(reshape(stnRec{ii,5}(:,2:13),1,[]))
        
    [cdfGrid, valGrid{ii}] = e_cdf(reshape(gridData{ii,5}(:,2:13),1,[]), 'bins', nBins);
    [cdfStn, valStn{ii}]   = e_cdf(reshape(  stnRec{ii,5}(:,2:13),1,[]), 'bins', nBins);
    
    iterMod = mod(ii-1, nSeries)+1;
    
    if iterMod == 1
        strSeries = cell(2*nSeries,1);
        hTs = nan(2*nSeries,1);
        
        hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'visible', 'off');
        hold on
    end
    
    if all(isnan(valGrid{ii})) || all(isnan(valStn{ii}))
       continue 
    end

        warning('off', 'fitness:unequalSize');
    stnBias(ii)  = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),  'bias');
    stnMAE(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'MAE');
    stnMAPE(ii)  = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),  'MAPE');
    stnWMAPE(ii) = fitness(valStn{ii}(2:end), valGrid{ii}(2:end), 'WMAPE');
    stnCorr(ii,:) = fitness(valStn{ii}(2:end), valGrid{ii}(2:end), 'pearson');
%     stnNse(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'nse');
%     stnKge(ii)   = fitness(valStn{ii}(2:end), valGrid{ii}(2:end),   'kge');
        warning('on', 'fitness:unequalSize');
    
    indCurr = [2*(iterMod-1) + 1, 2*iterMod];
    hTs(indCurr) = plot(...
         valStn{ii}(2:end),  cdfStn(2:end), '-', ...
         valGrid{ii}(2:end), cdfGrid(2:end), ':');
     
    %Set properties of current series:
    set(hTs(indCurr),'LineWidth',lnWd, ...
        'Color',colorsUse(iterMod,:));
     
    strSeries{indCurr(1)} = ['(' num2str(round2(stnRec{ii,3},2)) ', ' num2str(round2(stnRec{ii,2},2)) '): '  num2str(stnRec{ii,4})];
    strSeries{indCurr(2)} = ['(' num2str(round2(gridData{ii,3},2)) ', ' num2str(round2(gridData{ii,2},2)) ')'];
     
     
     if iterMod == nSeries || ii == length(gridData(:,1))  
        %Create legend:
        indRemLgd = find(isnan(hTs) == 1);
        hTs(indRemLgd) = [];
        strSeries(indRemLgd) = [];

        hLgd = legend(hTs, strSeries);
        
        hYLab = ylabel('Cumulative Probability');
        hXLab = xlabel([var ' (' units ')']);
        hTtl = title('Comparison of Station CDFs for GHCN Records Relative to Downscaled Data');
        
        
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
        cntr = cntr + 1;
        %Export:
        pathCdfPlot = fullfile(path, ...
            ['CDF_stn_', var, '_' num2str(cntr)]);
        %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        savefig(hFig,[pathCdfPlot '.fig']);
        print(hFig,[pathCdfPlot '.eps'],'-depsc2');
        print(hFig,[pathCdfPlot '.png'],'-dpng','-r600');
        
        
        hold off
     end
end
hold off

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

% 
% stats(5).name = 'Stn NSE';
% stats(5).data = stnNse;
% 
% stats(6).name = 'Stn KGE';
% stats(6).data = stnKge;

