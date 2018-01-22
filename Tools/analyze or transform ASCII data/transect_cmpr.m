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

%%SCRIPT TO COMPARE CROSSECTIONS OF ESRI FORMAT ASCII FILES
%Can handle both ESRI formatted ASCII files and WorldClim Binary file
%format.

%notes: 


clearvars -except
clc

%%Select Options:
metVar = 'Temperature (Deg C)';

lat  = 35.3; %Pakistan [60, 83, 21.5, 42] 30.5 is good
lonW = 70.1;
lonE = 79.5;

nCmpr = 2;  %Can compare 7 (the last is black!) at a time with current color vector ('ColorVec').
multAxis = 0;   %0 = NO OR 1 = yes; Series with Different Y axis must be the last entered.

colorVec = [0,0,1; 1,0,0; 0,1,1; 0,1,0; 1,0,1; 1,1,0; 0,0,0];





%%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
indParent = regexpi(pathScript, filesep);
indParent = indParent(end); %Assumes scripts are in distribution hierarchy.
addpath(genpath(pathScript(1:indParent-1))); %Add modules


%%INNERWORKINGS (DO NOT MODIFY):
if multAxis == 1
    disp('Elevation must be the second data product entered.');
end

%%User Prompt:
nameProduct = cell(nCmpr,1);
pathProduct = cell(nCmpr,1);
nScale      = cell(nCmpr,1);
lineStyle   = cell(nCmpr,1);

searchPath = pwd;
for ii = 1 : nCmpr
    uiwait(msgbox(sprintf('Select the climate product file to use. \n'), ... 
        '(Click OK to Proceed)','modal'));
    [fileTemp, pathTemp, ~]  = uigetfile({'*.asc;*.bil'},'Select the file to compute a transect of',searchPath);
        searchPath = pathTemp;
    pathProduct{ii} = fullfile(pathTemp,fileTemp);
    nameProduct{ii} = input('Enter the name of the climate product you selected:','s');
    lineStyle{ii} = input(['Enter the line style to use for this data series (' ...
        char(39) '-' char(39) ', ' char(39) '--' char(39) ', ' char(39) ':' ...
        char(39) ', or ' char(39) '-. ):'],'s');
    
    %Scale input data
    if regexpi(metVar,'pre')
       if ~isempty(regexpi(nameProduct{ii},'CRU')) && (~isempty(regexpi(nameProduct{ii},'0.5')) || ~isempty(regexpi(nameProduct{ii},'raw'))) %CRU pre is scaled by ten.
           nScale{ii} = 0.1;
       else     %W&M pre, PRISM pre, and WorldClim pre are not scaled
           nScale{ii} = 1;
       end
    elseif ~isempty(regexpi(metVar,'tmn')) || ~isempty(regexpi(metVar,'temp'))
       if ~isempty(regexpi(nameProduct{ii},'CRU')) && (~isempty(regexpi(nameProduct{ii},'0.5')) || ~isempty(regexpi(nameProduct{ii},'raw')))  %CRU temp is scaled by ten.
           nScale{ii} = 0.1;
       elseif ~isempty(regexpi(nameProduct{ii},'WorldClim')) || ~isempty(regexpi(nameProduct{ii},'WC')) %WorldClim temp is scaled by ten.
           nScale{ii} = 0.1;
       elseif regexpi(nameProduct{ii},'PRISM')  %PRISM temp is scaled by ten.
           nScale{ii} = 0.1;
       else     %W&M temp are not scaled
           nScale{ii} = 1;
       end  
    elseif ~isempty(regexpi(metVar,'unitless')) || ~isempty(regexpi(metVar,'anomaly'))
        nScale{ii} = 1;
    else
           error('Meteorological variable not recognized.');
    end
end





%%Process User Input Data:
coordinates = [lonW lonE lat];

tranOrg = cell(nCmpr,1);
transect = cell(nCmpr,1);
hdr      = cell(nCmpr,1);
interpFrom = cell(nCmpr,1);
orgRes = cell(nCmpr,1);
fineRes = 0;

%Find transect for each data file and define grid to use
for ii = 1 : nCmpr
[tranOrg{ii}, hdr{ii}] = climate_crossection(pathProduct{ii}, coordinates);
tranOrg{ii} = nScale{ii}*tranOrg{ii};
orgRes{ii} = hdr{ii}(1);

%Ensure 'tranOrg' has desired row 
yVec = ((hdr{ii}(2)-1)*hdr{ii}(5) + hdr{ii}(4) : -hdr{ii}(5) : hdr{ii}(4));
latInd = find(yVec <= lat);
latInd = latInd(1);
%Perform a check to ensure desired latitude is included in 'tranOrg':
if yVec(latInd) + hdr{ii}(5) >= lat && yVec(latInd) <= lat
   tranOrg{ii} = tranOrg{ii}(latInd,:);
else
    disp(['Error, latitude is not contained in tranOrg for data product ' num2str(ii) '.'])
end

    if orgRes{ii} > fineRes
        fineRes = orgRes{ii};
    end

    %Define grid using left-corner of box (Needs +1 to define x lower-right corner)
    interpFrom{ii} = linspace(hdr{ii}(3)+hdr{ii}(5)/2, hdr{ii}(3)+hdr{ii}(5)/2 + (orgRes{ii}-1)*hdr{ii}(5), orgRes{ii});
    interpFrom{ii} = [hdr{ii}(3), interpFrom{ii}, hdr{ii}(3)+hdr{ii}(5) + (orgRes{ii}-1)*hdr{ii}(5)];
    tranOrg{ii} = cat(2, tranOrg{ii}(:,1), tranOrg{ii}(:,:)  );
    tranOrg{ii} = cat(2, tranOrg{ii}(:,:), tranOrg{ii}(:,end));    
end

if fineRes < (lonE - lonW)*120  %If resolution is less than 30 arc-seconds, increase it.
   fineRes = (lonE - lonW)*120;
end

close all   %Closes any figures
hold off;   %Ensure hold is not on.
interpTo   = cell(nCmpr,1);

for ii = 1 : nCmpr
    interpTo{ii} = linspace(lonW, lonE, fineRes);
    
    transect{ii} = interp1(interpFrom{ii}, tranOrg{ii}(1,:), interpTo{ii} ,'nearest');
    
    if multAxis == 1 && ii == 1
        continue;
    elseif multAxis == 1 && ii == 2
        [AX,H1,H2] = plotyy(interpTo{1}, transect{1}, interpTo{2}, transect{2}, 'plot');

        set(get(AX(1),'Ylabel'),'String', metVar); 
        set(get(AX(2),'Ylabel'),'String','Elevation (m)');
        
        set(H1,'LineStyle', lineStyle{1}, 'LineWidth',5, 'Color',colorVec(1,:));
        set(H2,'LineStyle', lineStyle{2}, 'LineWidth',5, 'Color',colorVec(2,:));
        continue;
    else
        plot(interpTo{ii}, transect{ii}, lineStyle{ii}, 'LineWidth',5, ...
                                                    'Color',colorVec(ii,:));
    end
    
    if ii == 1
       hold on 
    end 
    %scatter(interpFrom{ii},tranOrg{ii}(1,:))
end



%{
%CHECK DELTA
check = (CRU_TS_cross ./ CRU_Norm_cross) .* WC_cross;
plot( (1:length(check(1,:))), check(1,:) ,'Color','black','LineWidth',2 );
%}

fSize = 40;
%%PLOT SPECS
%title(['Cross-Sectional Profile at Latitude ' num2str(lat) ' (degrees)']);
xlabel('Longitude (degrees)','FontSize',fSize,'FontWeight','bold');
ylabel(metVar,'FontSize',fSize,'FontWeight','bold');
hleg = legend(nameProduct,'Location','NorthEast','linewidth',4);
set(gca,'FontSize',fSize,'FontWeight','bold','linewidth',5);
%set(gca);
xlim([lonW lonE]);

hold off

%refline(0,0);

