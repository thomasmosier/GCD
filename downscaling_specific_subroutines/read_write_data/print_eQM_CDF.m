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

function print_eQM_CDF(sMP, sMH, sOH, sMeta, nBins, nPts, path)

error('printeQMCDF:obsolete','Update this function')
%Randomly choose grid cells for number of points specified:
%Doesn't pick cells along edges:


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

sMPBc = eQM_geodata(sMP, sMH, sOH);


sOH = ecdf_geodata(sOH, nBins, 1);
sMH = ecdf_geodata(sMH, nBins, 1);
sMP = ecdf_geodata(sMP, nBins, 1);
sMPBc = ecdf_geodata(sMPBc, nBins, 1);

%Remove points from consideration if data are nan at those locations:
iLat = NaN(nPts,1);
iLon = NaN(nPts,1);
nIter = 0;
for ii = 1 : nPts
    while (isnan(iLat(ii)) || isnan(iLon(ii))) && nIter < 200
        iLat(ii) = ceil(rand(1)*(length(sMH.data(1,:,1))-1));
        iLon(ii) = ceil(rand(1)*(length(sMH.data(1,1,:))-1));
        if numel(sOH.val{iLat(ii),iLon(ii)}) == 1 || numel(sMH.val{iLat(ii),iLon(ii)}) == 1 || numel(sMP.val{iLat(ii),iLon(ii)}) == 1
            iLat(ii) = NaN;
            iLon(ii) = NaN;
        end
        nIter = nIter + 1;
    end
end


indU = curr_ind(sMeta);

if regexpbl(sMeta.currVar,{'pr'})
    metVarDisp = 'Precipitation';
elseif regexpbl(sMeta.currVar,{'tmn'})
    metVarDisp = 'Minimum Temperature';
elseif regexpbl(sMeta.currVar,{'tmp'})
    metVarDisp = 'Mean Temperature';
elseif regexpbl(sMeta.currVar,{'tmx'})
    metVarDisp = 'Maximum Temperature';
else
    metVarDisp = 'Unknown';
end


%Make plot for each randomly generated point:
for ii = 1 : nPts
    hQMFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Visible','Off');
    
    hQMAxes = plot(...
         sOH.val{iLat(ii),iLon(ii)}(2:end),    sOH.cdf{iLat(ii),iLon(ii)}(2:end), ...
         sMH.val{iLat(ii),iLon(ii)}(2:end),    sMH.cdf{iLat(ii),iLon(ii)}(2:end), ...
        sMP.val{iLat(ii),iLon(ii)}(2:end),   sMP.cdf{iLat(ii),iLon(ii)}(2:end),  ...
      sMPBc.val{iLat(ii),iLon(ii)}(2:end), sMPBc.cdf{iLat(ii),iLon(ii)}(2:end));
    
    gcmValBC = (sOH.val{iLat(ii),iLon(ii)} - sMH.val{iLat(ii),iLon(ii)}) + sMP.val{iLat(ii),iLon(ii)} ;
    currRange = sMPBc.val{iLat(ii),iLon(ii)}(end)-sMPBc.val{iLat(ii),iLon(ii)}(1);
%     if nanmean(abs(gcmValBC - sMPBc.val{iLat(ii),iLon(ii)})) > 0.1*currRange
%         warning('print_downscale_cdf:bcDiff',['The bias-corrected data'...
%             ' appear to differ from the expected values by more than 10'...
%             ' percent of the data range for the pixel centered at ' num2str(iLat(ii)) 'N, ', num2str(iLon(ii)) 'E.'])
%     end
        
    ylim([sOH.cdf{iLat(ii),iLon(ii)}(2), sOH.cdf{iLat(ii),iLon(ii)}(end)]);
    set(gca,'FontSize',16, 'FontWeight','bold')

    hLgd = legend(hQMAxes(:), ...
                sMeta.hisObs{indU}, ...
                sMeta.hisMod{indU}, ...
                sMeta.projMod{indU}, ...
                [sMeta.projMod{indU} '(Bias-Corrected)']);
    hYLab = ylabel('Cumulative Probability');
    untBC = NC_units(sMP.attData);
    hXLab = xlabel([char(regexprep(metVarDisp,'(\<[a-z])','${upper($1)}')) ' (' ...
        char(untBC) ')']);
    hTtl = title(['CDFs for ' ...
        num2str(sMeta.currTime(2)) ' at Latitude = ' ...
        num2str(round2(sMP.lat(iLat(ii)),1)) ', Longitude = ' ...
        num2str(round2(sMP.lon(iLon(ii)),1))]);
    
    %Set figure data properties:
    set(hQMAxes(:), ...
        'LineWidth', lnWd);
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

    set(hQMAxes(1), ...
        'color', clrBrewer(1,:));
    set(hQMAxes(2), ...
        'color', clrBrewer(2,:));
    set(hQMAxes(3), ...
        'color', clrBrewer(3,:));
    set(hQMAxes(4), ...
        'color', clrBrewer(4,:));

    if ~exist(path,'dir')
        mkdir(path)
    end
    %Export:
    pathCdfPlot = fullfile(path, ...
        ['CDF_', sMeta.region, '_', num2str(iLat(ii)), '_', ...
        num2str(iLon(ii)), '_', sMeta.currVar, '_', ...
        num2str(sMeta.yrsOut(1)), '-', num2str(sMeta.yrsOut(2)), '_', ...
        num2str(sMeta.currTime(2))]);
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hQMFig,[pathCdfPlot '.fig']);
    print(hQMFig,[pathCdfPlot '.eps'],'-depsc2');
    print(hQMFig,[pathCdfPlot '.png'],'-dpng','-r600');
end

close(hQMFig);