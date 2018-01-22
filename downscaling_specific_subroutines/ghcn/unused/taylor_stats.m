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

function stats = taylor_stats(gridData, stnRec, nBins, metVar, dirReport)



for ii = 1 : length(gridData(:,1))
    S = allstats(reshape(stnRec{ii,5}(:,2:13),[],1),reshape(gridData{ii,5}(:,2:13),[],1));
    MYSTATS(ii,:) = S(:,2); % We get stats versus reference
end%for iserie
 MYSTATS(1,:) = S(:,1); % We assign reference stats to the first row

 keyboard
 szFrame = [12,7];
 
hFig = figure('Units','in','Position',[2 2 szFrame],'paperunits','in','paperposition',[2 2 szFrame]);
taylordiag(MYSTATS(1:5,2),MYSTATS(1:5,3),MYSTATS(1:5,4));




% [pp tt axl] = taylordiag(squeeze(statm(:,2)),squeeze(statm(:,3)),squeeze(statm(:,4)),...
%             'tickRMS',[25:25:150],'titleRMS',0,'tickRMSangle',135,'showlabelsRMS',0,'widthRMS',1,...
%             'tickSTD',[25:25:250],'limSTD',250,...
%             'tickCOR',[.1:.1:.9 .95 .99],'showlabelsCOR',1,'titleCOR',1);
% 
% for ii = 1 : length(tt)
%     set(tt(ii),'fontsize',9,'fontweight','bold')
%     set(pp(ii),'markersize',12)
%     if ii == 1
%         set(tt(ii),'String','Buoy');
%     else
%         set(tt(ii),'String',alphab(ii-1));
%     end
% end
% title(sprintf('%s: Taylor Diagram at CLIMODE Buoy','B'),'fontweight','bold');
% 
% tt = axl(2).handle;
% for ii = 1 : length(tt)
%     set(tt(ii),'fontsize',10,'fontweight','normal');
% end
% set(axl(1).handle,'fontweight','normal');



% 
% 
% for ii = 1 : length(gridData(:,1))
%     gridVec = cat(1, gridVec, reshape(gridData{ii,5}(:,2:13),[],1));
%     stnVec = cat(1, stnVec, reshape(stnData{ii,5}(:,2:13),[],1));
% end
% 
% [cdfGrid, valGrid] = ecdf_data(gridVec, nBins);
% [cdfStn, valStn]  = ecdf_data(stnVec, nBins);
% 
% aggBias  = fitness(valStn(2:end), valGrid(2:end),  'bias');
% aggMAE   = fitness(valStn(2:end), valGrid(2:end),   'MAE');
% aggMAPE  = fitness(valStn(2:end), valGrid(2:end),  'MAPE');
% aggWMAPE = fitness(valStn(2:end), valGrid(2:end), 'WMAPE');
% 
% 
% taylordiag(taylorStats);