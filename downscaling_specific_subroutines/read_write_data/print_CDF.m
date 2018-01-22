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

function print_CDF(path, sSimOut, sSim2Bc, sSimH, sRefH, var, yrsOut, yrsBase, nPts)

%Randomly choose grid cells for number of points specified:
%Doesn't pick cells along edges:

strVis = 'off';

varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';

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

% sRefH = ecdf_geodata(sRefH, nBins, 1);
% sSimH = ecdf_geodata(sSimH, nBins, 1);
% sSim2Bc = ecdf_geodata(sSim2Bc, nBins, 1);
% sSimOut = ecdf_geodata(sSimOut, nBins, 1);

%Set n points to be at most the number of points of data
nPts = min(nPts, numel(sSimOut.(varLat))*numel(sSimOut.(varLon)));

indLatSimOut = nan(nPts, 1);
indLonSimOut = nan(nPts, 1);
indLatSim2Bc = nan(nPts, 1);
indLonSim2Bc = nan(nPts, 1);
indLatSimH = nan(nPts, 1);
indLonSimH = nan(nPts, 1);
indLatRefH = nan(nPts, 1);
indLonRefH = nan(nPts, 1);

latPlt = nan(nPts, 1);
lonPlt = nan(nPts, 1);

nBrk = 10;

for ii = 1 : nPts
    cntrPt = 0;
    cntrBrk = 0;
    while cntrPt == 0
        indLatSimOut(ii) = randi(numel(sSimOut.(varLat)), 1);
        indLonSimOut(ii) = randi(numel(sSimOut.(varLon)), 1);
        
        latPlt(ii) = sSimOut.(varLat)(indLatSimOut(ii));
        lonPlt(ii) = sSimOut.(varLon)(indLonSimOut(ii));
        
        if any(~isnan(sSimOut.(var)(:,indLatSimOut(ii),indLonSimOut(ii))))
            cntrPt = 1;
            
            [~, indLatSim2Bc(ii)] = min(abs(sSim2Bc.(varLat) - latPlt(ii)));
            [~, indLonSim2Bc(ii)] = min(abs(sSim2Bc.(varLon) - lonPlt(ii)));
            
            [~, indLatSimH(ii)] = min(abs(sSimH.(varLat) - latPlt(ii)));
            [~, indLonSimH(ii)] = min(abs(sSimH.(varLon) - lonPlt(ii)));
            
            [~, indLatRefH(ii)] = min(abs(sRefH.(varLat) - latPlt(ii)));
            [~, indLonRefH(ii)] = min(abs(sRefH.(varLon) - lonPlt(ii)));
        end
        
        cntrBrk = cntrBrk + 1;
        
        if cntrBrk > nBrk
           break 
        end
    end
end



%Crop years:
%Find time-series elements to compare during bias correction:
indSimH = find(sSimH.(varDate)(:,1) >= min(yrsBase) & sSimH.(varDate)(:,1) <= max(yrsBase));
indRefH = find(sRefH.(varDate)(:,1) >= min(yrsBase) & sRefH.(varDate)(:,1) <= max(yrsBase));

if numel(indRefH) ~= numel(indSimH)
    [indSameSimH, ~] = ismember(sSimH.(varDate)(:,1), sRefH.(varDate)(:,1), 'rows');
    [indSameRefH, ~] = ismember(sRefH.(varDate)(:,1), sSimH.(varDate)(:,1), 'rows');
    
    indSimH = intersect(indSimH, sort(find(indSameSimH ~= 0)));
    indRefH = intersect(indRefH, sort(find(indSameRefH ~= 0)));
	warning('eQMgeodata:unequalYrs', ['The number of years present in the '...
       'two datasets being compared are not equal. They are going to be set to the common years of ' num2str(min(sSimH.(varDate)(indSimH,1))) '-' num2str(max(sSimH.(varDate)(indSimH,1)))]); 
end

%Find time-series elements to compare for output:
indSimOut = find(sSimOut.(varDate)(:,1) >= min(yrsOut) & sSimOut.(varDate)(:,1) <= max(yrsOut));
indSim2Bc = find(sSim2Bc.(varDate)(:,1) >= min(yrsOut) & sSim2Bc.(varDate)(:,1) <= max(yrsOut));

if numel(indSim2Bc) ~= numel(indSimOut)
    [indSameSimOut, ~] = ismember(sSimOut.(varDate)(:,1), sSim2Bc.(varDate)(:,1), 'rows');
    [indSameSim2Bc, ~] = ismember(sSim2Bc.(varDate)(:,1), sSimOut.(varDate)(:,1), 'rows');
    
    indSimOut = intersect(indSimOut, sort(find(indSameSimOut ~= 0)));
    indSim2Bc = intersect(indSim2Bc, sort(find(indSameSim2Bc ~= 0)));
	warning('eQMgeodata:unequalYrs', ['The number of years present in the '...
       'two datasets being compared are not equal. They are going to be set to the common years of ' num2str(min(sSimOut.(varDate)(indSimOut,1))) '-' num2str(max(sSimOut.(varDate)(indSimOut,1)))]); 
end

mnths = unique(sSimOut.(varDate)(:,2));
if numel(mnths) == 1
    strMnth = mnth_str(mnths);
elseif  numel(mnths) == 12
    strMnth = 'allMnths';
else
    strMnth = blanks(0);
    for ii = 1 : numel(mnths)
        strMnth = [strMnth, mnth_str(mnths(ii))];
    end
end


cdfSimOut = cell(nPts, 1);
cdfSim2Bc = cell(nPts, 1);
cdfRefH   = cell(nPts, 1);
cdfSimH   = cell(nPts, 1);
valSimOut = cell(nPts, 1);
valSim2Bc = cell(nPts, 1);
valRefH   = cell(nPts, 1);
valSimH   = cell(nPts, 1);

for ii = 1 : nPts
    [cdfSimOut{ii}, valSimOut{ii}] = e_cdf(sSimOut.(var)(indSimOut, indLatSimOut(ii), indLonSimOut(ii)));
    [cdfSim2Bc{ii}, valSim2Bc{ii}] = e_cdf(sSim2Bc.(var)(indSim2Bc, indLatSim2Bc(ii), indLonSim2Bc(ii)));
    [  cdfRefH{ii},   valRefH{ii}] = e_cdf(  sRefH.(var)(  indRefH,   indLatRefH(ii),   indLonRefH(ii)));
    [  cdfSimH{ii},   valSimH{ii}] = e_cdf(  sSimH.(var)(  indSimH,   indLatSimH(ii),   indLonSimH(ii)));
end
clear ii

%Remove points where insufficient data present
for ii = nPts : -1 : 1
    if all(isnan(valSimOut{ii})) || all(isnan(valSim2Bc{ii})) || all(isnan(valRefH{ii})) || all(isnan(valSimH{ii}))
        cdfSimOut(ii) = [];
        cdfSim2Bc(ii) = [];
        cdfRefH(ii) = [];
        cdfSimH(ii) = [];
        valSimOut(ii) = [];
        valSim2Bc(ii) = [];
        valRefH(ii) = [];
        valSimH(ii) = [];
        
        indLatSimOut(ii) = [];
        indLonSimOut(ii) = [];
        indLatSim2Bc(ii) = [];
        indLonSim2Bc(ii) = [];
        indLatSimH(ii) = [];
        indLonSimH(ii) = [];
        indLatRefH(ii) = [];
        indLonRefH(ii) = [];

        latPlt(ii) = [];
        lonPlt(ii) = [];
    end
end
clear ii

nPts = numel(cdfSimOut(:));

%Return if no points remaining
if nPts == 0
   return 
end


if regexpbl(var,{'pr'})
    metVarDisp = 'Precipitation';
    unit = 'mm';
elseif regexpbl(var,{'tmn'})
    metVarDisp = 'Minimum Temperature';
    unit = '\circ C';
elseif regexpbl(var,{'tmp'})
    metVarDisp = 'Mean Temperature';
    unit = '\circ C';
elseif regexpbl(var,{'tmx'})
    metVarDisp = 'Maximum Temperature';
    unit = '\circ C';
else
    metVarDisp = 'Unknown';
    unit = '?';
end


%Make plot for each randomly generated point:
for ii = 1 : nPts
    hQMFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame], 'Visible', strVis);
    
    hQMAxes = plot(...
          valSimH{ii},   cdfSimH{ii}, ...
      	  valRefH{ii},   cdfRefH{ii}, ...
        valSim2Bc{ii}, cdfSim2Bc{ii}, ...
        valSimOut{ii}, cdfSimOut{ii});
    
    yRang = [min([min(cdfSimH{ii}), min(cdfRefH{ii}), min(cdfSim2Bc{ii}), min(cdfSimOut{ii})]), max([max(cdfSimH{ii}), max(cdfRefH{ii}), max(cdfSim2Bc{ii}), max(cdfSimOut{ii})])];
        
    ylim(yRang);
%     set(gca,'FontSize',16)

    hLgd = legend(hQMAxes(:), ...
                'Simulation (Historical)', ...
                'Reference (Historical)', ...
                'Simulation (Projection, Input)', ...
                'Simulation (Projection, Output)', ...
                'location', 'southeast');
    hYLab = ylabel('Cumulative Probability (%)');
    hXLab = xlabel([char(regexprep(metVarDisp,'(\<[a-z])','${upper($1)}')) ' (' ...
        char(unit) ')']);
    hTtl = title(['CDFs for ' ...
        strMnth ' at Latitude = ' ...
        num2str(round2(sSimOut.(varLat)(indLatSimOut(ii)),1)) ', Longitude = ' ...
        num2str(round2(sSimOut.(varLon)(indLonSimOut(ii)),1))]);
    
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


    %WRITE TO FILE
    if ~exist(path,'dir')
        mkdir(path)
    end
    %Export:
    pathCdfPlot = fullfile(path, ...
        ['cdf_', num2str(indLatSimOut(ii)), '_', ...
        num2str(indLatSimOut(ii)), '_', var, ...
        num2str(yrsOut(1)), '-', num2str(yrsOut(2)), '_', ...
        strMnth]);
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hQMFig,[pathCdfPlot '.fig']);
    print(hQMFig,[pathCdfPlot '.eps'],'-depsc2');
    print(hQMFig,[pathCdfPlot '.png'],'-dpng','-r600');
end

close(hQMFig);