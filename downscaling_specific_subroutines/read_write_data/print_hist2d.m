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

function print_hist2d(path, sSimOut, sSimOutV, sSim2Bc, sSim2BcV, sSimH, sSimHV, sRefH, sRefHV, var, varV, yrsOut, yrsBase, nPts)


%Randomly choose grid cells for number of points specified:
%Doesn't pick cells along edges:

nBinsH = 20;
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
szFrame = [7,5];


%Initialize arrays:
pearSimOut = nan(nPts, 1);
pearSim2Bc = nan(nPts, 1);
pearRefH = nan(nPts, 1);
pearSimH = nan(nPts, 1);

nSimOut = cell(nPts,1);
nSim2Bc = cell(nPts,1);
nRefH = cell(nPts,1);
nSimH = cell(nPts,1);

cntrs = cell(nPts,1);
edgs1 = cell(nPts,1);
edgs2 = cell(nPts,1);

% ndMPBC = nan(numel(sRefHV.(varLat)),numel(sRefHV.(varLon)));
% ndMP = nan(numel(sRefHV.(varLat)),numel(sRefHV.(varLon)));
% ndMH = nan(numel(sRefHV.(varLat)),numel(sRefHV.(varLon)));

indlatSimOut = nan(nPts, 1);
indlonSimOut = nan(nPts, 1);
indLatSim2Bc = nan(nPts, 1);
indLonSim2Bc = nan(nPts, 1);
indLatSimH = nan(nPts, 1);
indLonSimH = nan(nPts, 1);
indLatRefH = nan(nPts, 1);
indLonRefH = nan(nPts, 1);

latPlt = nan(nPts, 1);
lonPlt = nan(nPts, 1);

for ii = 1 : nPts
    cntrPt = 0;
    while cntrPt == 0
        indlatSimOut(ii) = randi(numel(sSimOut.(varLat)), 1);
        indlonSimOut(ii) = randi(numel(sSimOut.(varLon)), 1);
        
        latPlt(ii) = sSimOut.(varLat)(indlatSimOut(ii));
        lonPlt(ii) = sSimOut.(varLon)(indlonSimOut(ii));
        
        if any(~isnan(sSimOut.(var)(:,indlatSimOut(ii),indlonSimOut(ii))))
            cntrPt = 1;
            
            [~, indLatSim2Bc(ii)] = min(abs(sSim2Bc.(varLat) - latPlt(ii)));
            [~, indLonSim2Bc(ii)] = min(abs(sSim2Bc.(varLon) - lonPlt(ii)));
            
            [~, indLatSimH(ii)] = min(abs(sSimH.(varLat) - latPlt(ii)));
            [~, indLonSimH(ii)] = min(abs(sSimH.(varLon) - lonPlt(ii)));
            
            [~, indLatRefH(ii)] = min(abs(sRefH.(varLat) - latPlt(ii)));
            [~, indLonRefH(ii)] = min(abs(sRefH.(varLon) - lonPlt(ii)));
        end
    end
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


for ii = 1 : nPts
    matSimPOut = [reshape(sSimOut.(var)(indSimOut,indlatSimOut(ii),indlonSimOut(ii)),[],1), reshape(sSimOutV.(varV)(indSimOut,indlatSimOut(ii),indlonSimOut(ii)),[],1)];
    matSimP2Bc = [reshape(sSim2Bc.(var)(indSim2Bc,indLatSim2Bc(ii),indLonSim2Bc(ii)),[],1), reshape(sSim2BcV.(varV)(indSim2Bc,indLatSim2Bc(ii),indLonSim2Bc(ii)),[],1)];
    matSimH    = [reshape(  sSimH.(var)(  indSimH,  indLatSimH(ii),  indLonSimH(ii)),[],1), reshape(  sSimHV.(varV)(  indSimH,  indLatSimH(ii),  indLonSimH(ii)),[],1)];
    matRefH    = [reshape(  sRefH.(var)(  indRefH,  indLatRefH(ii),  indLonRefH(ii)),[],1), reshape(  sRefHV.(varV)(  indRefH,  indLatRefH(ii),  indLonRefH(ii)),[],1)];
    
    if all(all(~isnan(matSimPOut)))
        %Ensure temperature first:
        if regexpbl(varV,'pr')
            matSimPOut = matSimPOut(:,[2,1]);
            matSimP2Bc = matSimP2Bc(:,[2,1]);
            matSimH = matSimH(:,[2,1]);
            matRefH = matRefH(:,[2,1]);
        end

        rng1(1) = min([matSimPOut(:,1);matSimP2Bc(:,1);matSimH(:,1);matRefH(:,1)]);
        rng2(1) = min([matSimPOut(:,2);matSimP2Bc(:,2);matSimH(:,2);matRefH(:,2)]);
        rng1(2) = max([matSimPOut(:,1);matSimP2Bc(:,1);matSimH(:,1);matRefH(:,1)]);
        rng2(2) = max([matSimPOut(:,2);matSimP2Bc(:,2);matSimH(:,2);matRefH(:,2)]);

        edgs1{ii} = linspace(rng1(1),rng1(2),nBinsH+1);
        edgs2{ii} = linspace(rng2(1),rng2(2),nBinsH+1);

        cntrs{ii} = {(0.5*diff(edgs1{ii}) + edgs1{ii}(1:end-1))',(0.5*diff(edgs2{ii}) + edgs2{ii}(1:end-1))'};

        nSimOut{ii} = hist3(matSimPOut,cntrs{ii}); 
          nSim2Bc{ii} = hist3(  matSimP2Bc,cntrs{ii}); 
          nRefH{ii} = hist3(  matRefH,cntrs{ii}); 
          nSimH{ii} = hist3(  matSimH,cntrs{ii});

        tauMPBC = corr(matSimPOut);
            pearSimOut(ii) = tauMPBC(2);
        tauMP = corr(matSimP2Bc);
            pearSim2Bc(ii) = tauMP(2);
        tauMH = corr(matSimH);
            pearSimH(ii) = tauMH(2);
        tauOH = corr(matRefH);
            pearRefH(ii) = tauOH(2);
    end 
end


tickPts = round(linspace(1,nBinsH,10));

%Make plot for each randomly generated point:
for ii = 1 : nPts
    %Make four subplots:
    hQMFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame],'visible', strVis);

    try %Put inside 'try' in case point selection did not work as intended and values are nan.
        %Historic Simulation:
        hSub1 = subplot(2,2,1);
        pcolor((1:nBinsH), (1:nBinsH), nSimH{ii});shading flat
        hXLab{1} = xlabel('Temperature (Deg C)'); hYLab{1} = ylabel('Precipitation (mm)'); 
        hTtl{1} = title(['Sim Hist (r = ' num2str(round2(pearSimH(ii),2)) ')']);
        ax = gca;
        ax.XTick = tickPts;
        ax.YTick = tickPts;
        ax.XTickLabel = round2(cntrs{ii}{1}(tickPts),0);
        ax.YTickLabel = round2(cntrs{ii}{2}(tickPts),0);


        %Historic Reference:
        hSub2 = subplot(2,2,2);
        pcolor((1:nBinsH), (1:nBinsH), nRefH{ii});shading flat
        hXLab{2} = xlabel('Temperature (Deg C)'); hYLab{2} = ylabel('Precipitation (mm)'); 
        hTtl{2} = title(['Ref Hist (r = ' num2str(round2(pearRefH(ii),2)) ')']);
        ax = gca;
        ax.XTick = tickPts;
        ax.YTick = tickPts;
        ax.XTickLabel = round2(cntrs{ii}{1}(tickPts),0);
        ax.YTickLabel = round2(cntrs{ii}{2}(tickPts),0);


        %Input Projection Simulation:
        hSub3 = subplot(2,2,3);
        pcolor((1:nBinsH), (1:nBinsH), nSim2Bc{ii});shading flat
        hXLab{3} = xlabel('Temperature (Deg C)'); hYLab{3} = ylabel('Precipitation (mm)'); 
        hTtl{3} = title(['Sim Proj (In) (r = ' num2str(round2(pearSim2Bc(ii),2)) ')']);
        ax = gca;
        ax.XTick = tickPts;
        ax.YTick = tickPts;
        ax.XTickLabel = round2(cntrs{ii}{1}(tickPts),0);
        ax.YTickLabel = round2(cntrs{ii}{2}(tickPts),0);


        %Output Projection Simulation:
        hSub4 = subplot(2,2,4);
        pcolor((1:nBinsH), (1:nBinsH), nSimOut{ii});shading flat
        hXLab{4} = xlabel('Temperature (Deg C)'); hYLab{4} = ylabel('Precipitation (mm)');
        hTtl{4} = title(['Sim Proj (Out) (r = ' num2str(round2(pearSimOut(ii),2)) ')']);
        ax = gca;
        ax.XTick = tickPts;
        ax.YTick = tickPts;
        ax.XTickLabel = round2(cntrs{ii}{1}(tickPts),0);
        ax.YTickLabel = round2(cntrs{ii}{2}(tickPts),0);
    
    catch ME
        if (strcmp(ME.identifier,'MATLAB:pcolor:InputSizeMismatch'))
            continue
        end
    end
    
    %Edit colorbar
    hCbAx = axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
    maxP =  max([max2d(nSimOut{ii}),max2d(nSim2Bc{ii}),max2d(nSimH{ii}),max2d(nRefH{ii})]);
    cbT = round2(linspace(0, maxP, min([maxP+1,8])),0);
    hCb = colorbar ('location','southoutside','FontSize',ftSz,'Ticks',cbT,'TickLabels',cbT);
    caxis(hSub1, [0, maxP]);
    caxis(hSub2, [0, maxP]);
    caxis(hSub3, [0, maxP]);
    caxis(hSub4, [0, maxP]);
    caxis(hCbAx, [0, maxP]);

    
    ylabel(hCb,'Joint Probability Density (%)');
    cpos = hCb.Position;
    cpos(4) = 0.6*cpos(4);
    cpos(2) = cpos(2)-0.01; 
    hCb.Position = cpos;
    
    %Change subfig positon:
    s1Pos = get(hSub1,'position'); % [left bottom width height]
    s2Pos = get(hSub2,'position'); % [left bottom width height]
    s3Pos = get(hSub3,'position'); % [left bottom width height]
    s4Pos = get(hSub4,'position'); % [left bottom width height]
    
    nudge = 0.05;
    s1Pos(3) = s1Pos(3) - nudge;
    s2Pos(3) = s2Pos(3) - nudge;
    s3Pos(3) = s3Pos(3) - nudge;
    s4Pos(3) = s4Pos(3) - nudge;
    s1Pos(4) = s1Pos(3);
    s2Pos(4) = s1Pos(3);
    s3Pos(4) = s1Pos(3);
    s4Pos(4) = s1Pos(3);
    s1Pos(2) = s1Pos(2) + 1.8*nudge;
    s2Pos(2) = s2Pos(2) + 1.8*nudge;
    s3Pos(2) = s3Pos(2) + 2.5*nudge;
    s4Pos(2) = s4Pos(2) + 2.5*nudge;
    
    %Old version
%     s1Pos(4) = s1Pos(4) - s1Pos(3);
%     s1Pos(2) = s1Pos(2) +( s1Pos(4) - s1Pos(3)) + 0.1;
% 
%     s2Pos(2) = s1Pos(2);
%     s2Pos(4) = s1Pos(4);
%     s2Pos(1) = s2Pos(1) + 0.05;

    set(hSub1,'position',s1Pos);
    set(hSub2,'position',s2Pos);
    set(hSub3,'position',s3Pos);
    set(hSub4,'position',s4Pos);
    
    
    
    %Edit Axes Properties:
    set([hXLab{:}, hYLab{:}]  , ...
        'FontSize'   , ftSzAx, ...
        'Color', 'black', ...
        'FontName'   , strFont);
    set([hTtl{:}], ...
        'FontWeight' , 'bold', ...
        'FontSize'   , ftSzAx);
    set([hSub1, hSub2, hSub3, hSub4], ...
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

    
    
    %WRITE TO FILE
    if ~exist(path,'dir')
        mkdir(path)
    end
    %Export:
    pathCdfPlot = fullfile(path, ...
        ['jointhist_', num2str(indlatSimOut(ii)), '_', ...
        num2str(indlatSimOut(ii)), '_', var, 'by', varV, ...
        num2str(yrsOut(1)), '-', num2str(yrsOut(2)), '_', ...
        strMnth]);
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    savefig(hQMFig,[pathCdfPlot '.fig']);
    print(hQMFig,[pathCdfPlot '.eps'],'-depsc2');
    print(hQMFig,[pathCdfPlot '.png'],'-dpng','-r600');
end

close(hQMFig);