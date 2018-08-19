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

function GHCN_2_grid_cmpr_v3(dirTs, sMeta, ghcnType, nScale, varargin)



%This script uses gridded data with ESRI header format and output from
%'GHCN_find', which are compared using r-sqaured and RMSE metrics.

%varargin{1} = GHCN ID of any stations to exclude from analysis

%Notes: 
%This script presumes that all ASCII files in the folder have the
%same naming scheme.  The only naming requirement is that the gridded 
%time-series files have a name of the form "*_year_month.asc".


yrs = sMeta.yrsOut;

if isempty(nScale)
    nScale = 1;
elseif nScale ~= 1
    warning('GHCN_2_grid_cmpr:scale',['The data is going to be scaled by a factor of ' num2str(nScale) '.  The scaling factor was an input parameter.'])
end

%%FIND GHCN STATIONS IN REGION OF INTEREST AND FIND GRIDDED DATA

%Use later:
%uiwait(msgbox(sprintf('Select the WorldClim Elevation file in binary format. \n'), ... 
%        '(Click OK to Proceed)','modal'));


if iscell(dirTs)
   nVar = numel(dirTs);
else
    nVar = 1;
end

metVar = '';

stnRec     = cell(nVar,1);
gridReform = cell(nVar,1);

if ~regexpbl(class(dirTs),'cell')
    temp = dirTs;
    dirTs = cell(1);
    dirTs{1} = temp;
end

%Attempt to put temperature first:
if numel(dirTs) == 2
   if regexpbl(dirTs{1},'pre') && regexpbl(dirTs{2},'tmp')
      dirTs = [dirTs(2); dirTs(1)] ;
   end
end
stnYrs = cell(nVar,1);

%LOOP OVER VARIABLES:
for kk = 1 : nVar
    if regexpbl(dirTs{kk},'tmp')
        metVar = 'tmp';
        units = 'celsius';
        sMeta.currVar = metVar;
    elseif regexpbl(dirTs{kk},'tmx')
        metVar = 'tmx';
        units = 'celsius';
        sMeta.currVar = metVar;
    elseif regexpbl(dirTs{kk},'tmn')
        metVar = 'tmn';
        units = 'celsius';
        sMeta.currVar = metVar;
    elseif regexpbl(dirTs{kk},'pre')
        units = 'mm';
        metVar = 'pre';
        sMeta.currVar = metVar;
    else
        continue
    end

    %Remove crd field and read_geodata will simply read all data present
    %(wont select bounding box)
    if isfield(sMeta,'crd')
        sMeta = rmfield(sMeta,'crd');
    end
    sGridded = read_geodata_v2(dirTs{kk}, metVar, nan(1,2), nan(1,2), nan(1,2), nan, nan);
    %sGridded = read_geodata(dirTs{kk}, sMeta, 'ts');
    %date = days_2_date(sGridded.time,[1900,1,1],'gregorian')
    %squeeze(sGridded.data(2,:,:))
    
    latTemp = box_edg(sGridded.latitude);
    lonTemp  = box_edg(sGridded.longitude);

    if regexpbl(ghcnType,'joint') 
        if regexpbl(metVar,'tmp')
            ghcnCurr = 'adj';
        else
            ghcnCurr = 'non';
        end
    else
        ghcnCurr = ghcnType;
    end
    %Find corresponding GHCN data:
    stnRec{kk} = GHCN_find([lonTemp(1), lonTemp(end), latTemp(end), latTemp(1)], yrs, metVar, ghcnCurr);
        %stnData is a cell array where: 
        %{i,1} is the GHCN station ID,
        %{i,2} is the station latitude, 
        %{i,3} is the station longitude,
        %{i,4} is the station elevation, and 
        %{i,5} is an m by 13 cell-array with yearly data (first column is year).

    %Inform user of level of GHCN data availability:
    if ~isempty(regexpi(ghcnType,'adj'))
        dataConMessage = 'cleaned';
    elseif regexpbl(ghcnType,{'org','orig','non'})
        dataConMessage = 'raw';
    elseif regexpbl(ghcnType,'joint')
        dataConMessage = 'raw precipitation and cleaned temperature';
    else
        error('GHCN_Grid_Cmpr:dataType',['The GHCN data string is ' ...
            char(39) ghcnType char(39) ', which is an unkown identifier.']);
    end    

    if isempty(stnRec{kk}(:,1))    %If no stn data, exit this function call.
        warning('GHCN_2_grid_cmpr:noStn',['Stats for ' strrep(dirTs{kk},':',filesep) ' using ' dataConMessage ' GHCN data are '...
            'not being computed because there are no stations in this region.']);
       return; 
    end

    %Find the number of years of data at each GHCN station:
    stnYrs{kk} = zeros(length(stnRec{kk}(:,1)),1);
    for ii = 1 : length(stnRec{kk}(:,1))
        stnYrs{kk}(ii) = length(stnRec{kk}{1,5}(:,1));
    end




    %%FIND GRIDDED DATA AT GHCN STATION LOCATIONS:
    %initialize cell-array with same dimensions as GHCN array, but NaN values:

    %  [dateRef, gcmUnits] = NC_time_units(attTime);
    %  datesGrid = days_2_date(sGridded.time, dateRef, NC_cal(attTime));

    gridReform{kk} = init_grid_GHCN(stnRec{kk});
    gridReform{kk} = geodata_2_GHCN_form(sGridded, gridReform{kk}, stnRec{kk}, 'var', metVar);

    % FOR DIAGNOSTICS:
    % for ii = 1 : numel(gridReform{kk}(:,1))
    %     numel(find(~isnan(gridReform{kk}{ii,5}(:,2:13))))
    % %     numel(find(~isnan(stnRec{kk}{ii,5}(:,2:13))))
    % end


    %Exclude GHCN stations if listed in preamble of script.
    if ~isempty(varargin) && ~isempty(varargin{1})
        [stnEx, gridReform{kk}, stnRec{kk}, stnYrs{kk}] = stn_exclude(gridReform{kk}, stnRec{kk}, stnYrs{kk}, varargin{1});
    else
        stnEx = [];
    end

    %Exclude GHCN stations if all gridded data is NaN at that location:
    [recStnNan, gridReform{kk}, stnRec{kk}, stnYrs{kk}] = del_ghcn_nan(gridReform{kk}, stnRec{kk}, stnYrs{kk});

    disp([ num2str(length(stnRec{kk}(:,1))) ' gauge sites with ' ...
        dataConMessage ' GHCN data are available for the period' ... 
        ' of interest within the geogrpahic region specified, with a total' ...
        ' of ' num2str(sum(stnYrs{kk})) ' station-years of data.']);

    %Check for negative precip values:
    if regexpi(metVar, 'pre')   %Check if it is precipitation and there are negative values.  If so, return error:
        for ii = 1 : length(gridReform{kk}(:,1))
            if ~isempty(find(gridReform{kk}{ii,5}(:,2:13) < 0, 1))
                warning('GHCN_2_grid_cmpr:negPreGrid',['Precipitation data in ' dirTs{kk} ' contains negative values.']);
            elseif ~isempty(find(stnRec{kk}{ii,5}(:,2:13) < 0, 1)) 
                warning('GHCN_2_grid_cmpr:negPreGHCN',['GHCN precipitation data for region associated with ' dirTs{kk} ' contains negative values.']);
            end
        end
    end

end

%No metVariable loaded. Don't go any further.
if isempty(metVar)
    warning('GHCN_2_grid_cmpr:wrongVar',['GHCN analysis is terminating '...
        'early because the current data are not mean temperature or precipitation.']);
    return
end

if regexpbl(ghcnType,'joint')
   nVar = 1; 
end
%%CALCULATE STATISTICS
for kk = 1 : nVar
    if regexpbl(ghcnType,'joint')
        %Create folder to write to:
        indSep = regexpi(dirTs{1},filesep);
        if regexpbl(dirTs{1},'tmp')
            indDir = 1;
        elseif regexpbl(dirTs{2},'tmp')
            indDir = 2;
        else
            indDir = 1;
            warning('GHCN_2_grid_cmpr:noTmp',['Mean Temperature is '...
                'expected as an input for joint variable analysis, but '...
                'it appears to not be present.']);
        end

        dirReport = fullfile(dirTs{indDir}(1:indSep(end-1)-1), ...
            ['GHCN_joint_analysis_', strrep(dirTs{indDir}(indSep(end-1)+1:indSep(end)-1),'_tmp','')]);
        if ~exist(dirReport,'dir')
            mkdir(dirReport)
        end
        disp(['Joint GHCN stats are being written to ' char(39) dirReport char(39)]);

        %Keep only stations that report both precip and temp:
        stnCrd{1} = [cell2mat(squeeze(stnRec{1}(:,2))),  cell2mat(squeeze(stnRec{1}(:,3)))];
        stnCrd{2} = [cell2mat(squeeze(stnRec{2}(:,2))),  cell2mat(squeeze(stnRec{2}(:,3)))];
        
        if isempty(stnCrd{1}) || isempty(stnCrd{2})
            warning('GHCN_2_grid_cmpr:noStns','Analysis is aborting because zero stations were found for one of the variables.');
            return
        end

        gridCrd{1} = [cell2mat(squeeze(gridReform{1}(:,2))),  cell2mat(squeeze(gridReform{1}(:,3)))];
        gridCrd{2} = [cell2mat(squeeze(gridReform{2}(:,2))),  cell2mat(squeeze(gridReform{2}(:,3)))];


    %     [stnUseID, ind1StnUse, ind2StnUse] = intersect(stnID{1},stnID{2});
        [~, ind1StnUse, ind2StnUse] = intersect(round2(stnCrd{1},1),round2(stnCrd{2},1),'rows');
        [~, ind1GridUse, ind2GridUse] = intersect(round2(gridCrd{1},1),round2(gridCrd{2},1),'rows');
        
        if ~isequal(ind1StnUse,ind1GridUse) || ~isequal(ind2StnUse,ind2GridUse) 
%             [~, ind1GridTest, ind1StnTest] = intersect(round2(gridCrd{1}(ind1GridUse,:),1),round2(stnCrd{1}(ind1StnUse,:),1),'rows');
%             ind1GridUse = ind1GridUse(ind1GridTest);
%             ind1StnUse = ind1StnUse(ind1StnTest);
            
            [~, ind1GridSame, ind1StnSame] = intersect(round2(gridCrd{1},1),round2(stnCrd{1},1),'rows');
            [~, ind2GridSame, ind2StnSame] = intersect(round2(gridCrd{2},1),round2(stnCrd{2},1),'rows');
            stnRec{1} =     stnRec{1}(ind1StnSame,:);
            stnRec{2} =     stnRec{2}(ind2StnSame,:);
            gridReform{1} =     gridReform{1}(ind1GridSame,:);
            gridReform{2} =     gridReform{2}(ind2GridSame,:);
            gridCrd{1} = gridCrd{1}(ind1GridSame,:);
            gridCrd{2} = gridCrd{2}(ind2GridSame,:);
            stnCrd{1} = stnCrd{1}(ind1StnSame,:);
            stnCrd{2} = stnCrd{2}(ind2StnSame,:);
            
            [~,  ind1StnUse,  ind2StnUse] = intersect(round2( stnCrd{1},1), round2( stnCrd{2},1),'rows');
            [~, ind1GridUse, ind2GridUse] = intersect(round2(gridCrd{1},1), round2(gridCrd{2},1),'rows');
            
            if ~isequal(ind1StnUse,ind1GridUse) || ~isequal(ind2StnUse,ind2GridUse) 
                warning('GHCN_2_grid_cmpr:jointStns','The station and grid indices for joint variable comparison do not line up. This could likely be solved through additional coding.') 
                return
            end
        end

            stnRec{1} = stnRec{1}(ind1StnUse,:);
            stnRec{2} = stnRec{2}(ind2StnUse,:);
            stnYrs{1} = stnYrs{1}(ind1StnUse);
            stnYrs{2} = stnYrs{2}(ind2StnUse);
        gridReform{1} = gridReform{1}(ind1StnUse,:);
        gridReform{2} = gridReform{2}(ind2StnUse,:);

        %Keep only time-series elements present for both precip and temp:
        for ii = numel(stnRec{1}(:,1)) : -1 : 1
            yrsStn{1} = stnRec{1}{ii,5}(:,1);
            yrsStn{2} = stnRec{2}{ii,5}(:,1);

            yrsAll = union(yrsStn{1},yrsStn{2});
            
            for jj = numel(yrsAll) : -1 : 1
                if isempty(find(yrsAll(jj) == stnRec{2}{ii,5}(:,1) ,1))
                        stnRec{1}{ii,5}(    stnRec{1}{ii,5}(:,1) == yrsAll(jj),:) = [];
                    gridReform{1}{ii,5}(gridReform{1}{ii,5}(:,1) == yrsAll(jj),:) = [];
                    stnYrs{1}(ii) = stnYrs{1}(ii) - 1;
                    continue
                elseif isempty(find(yrsAll(jj) == stnRec{1}{ii,5}(:,1) ,1))
                        stnRec{2}{ii,5}(    stnRec{2}{ii,5}(:,1) == yrsAll(jj),:) = [];
                    gridReform{2}{ii,5}(gridReform{2}{ii,5}(:,1) == yrsAll(jj),:) = [];
                    stnYrs{2}(ii) = stnYrs{2}(ii) - 1;
                    continue
                end

                indNan1 = find(isnan( stnRec{1}{ii,5}(stnRec{1}{ii,5}(:,1) == yrsAll(jj),:) ));
                indNan2 = find(isnan( stnRec{2}{ii,5}(stnRec{2}{ii,5}(:,1) == yrsAll(jj),:) ));

                if ~isempty(indNan1)
                        stnRec{2}{ii,5}(    stnRec{2}{ii,5}(:,1) == yrsAll(jj),indNan1) = nan;
                    gridReform{2}{ii,5}(gridReform{2}{ii,5}(:,1) == yrsAll(jj),indNan1) = nan;
                end
                if ~isempty(indNan2)
                        stnRec{1}{ii,5}(    stnRec{1}{ii,5}(:,1) == yrsAll(jj),indNan2) = nan;
                    gridReform{1}{ii,5}(gridReform{1}{ii,5}(:,1) == yrsAll(jj),indNan2) = nan;
                end
            end
        end
        

        %Adjust years used in each variable and between GHCN and gridded
        %data:
        for ii = numel(stnRec{1}(:,1)) : -1 : 1
            yrs1Stn = stnRec{1}{ii,5}(:,1);
            yrs1Grid = gridReform{1}{ii,5}(:,1);
            if ~isequal(yrs1Stn,yrs1Grid)
                [~, ind1Stn, ind1Grid] = intersect(yrs1Stn,yrs1Grid);
                if ~isempty(ind1Stn)
                    stnRec{1}{ii,5} = stnRec{1}{ii,5}(ind1Stn,:);
                    gridReform{1}{ii,5} = gridReform{1}{ii,5}(ind1Grid,:);
                else
                    stnRec{1}(ii,5) = [];
                    gridReform{1}(ii,5) = [];
                end
            end
            yrs2Stn = stnRec{2}{ii,5}(:,1);
            yrs2Grid = gridReform{2}{ii,5}(:,1);
            if ~isequal(yrs2Stn,yrs2Grid)
                [~, ind2Stn, ind2Grid] = intersect(yrs2Stn,yrs2Grid);
                if ~isempty(ind2Stn)
                    stnRec{2}{ii,5} = stnRec{1}{ii,5}(ind2Stn,:);
                    gridReform{2}{ii,5} = gridReform{2}{ii,5}(ind2Grid,:);
                else
                    stnRec{2}(ii,50) = [];
                    gridReform{2}(ii,5) = [];
                end
            end

            yrs1Stn = stnRec{1}{ii,5}(:,1);
            yrs2Stn = stnRec{2}{ii,5}(:,1);
            [~, ind1Stn, ind2Stn] = intersect(yrs1Stn,yrs2Stn);
            if ~isempty(ind1Stn)
                stnRec{1}{ii,5} = stnRec{1}{ii,5}(ind1Stn,:);
                stnRec{2}{ii,5} = stnRec{2}{ii,5}(ind2Stn,:);
                gridReform{1}{ii,5} = gridReform{1}{ii,5}(ind1Stn,:);
                gridReform{2}{ii,5} = gridReform{2}{ii,5}(ind2Stn,:);
            else
                stnRec{1}(ii,5) = [];
                stnRec{2}(ii,5) = [];
                gridReform{1}(ii,5) = [];
                gridReform{2}(ii,5) = [];
            end
        end
        
        nBins = 20;

    
        if numel(stnRec{1}(:)) == 0 || numel(stnRec{2}(:)) == 0 || numel(gridReform{1}(:)) == 0 || numel(gridReform{2}(:)) == 0
            warning('GHCN_2_grid_cmpr:dimmismatch','Joint GHCN analysis is being aborted because data dimenions are not the same.');
            return
        end
        
        [aggStats,   nCorrAggCom]  = agg_joint_stats(gridReform, stnRec, nBins, dirReport);
        [mnthStats, nCorrMnthCom] = mnth_joint_stats(gridReform, stnRec, nBins, dirReport);
        [stnStats,   nCorrStnCom]  = stn_joint_stats(gridReform, stnRec, nBins, dirReport);
        
        strType = 'joint histograms';
        fileGeoStats = ['GHCN_joint_hist_stn_geoStats_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.asc'];
        nameReport = ['GHCN_joint_hist_tableStats_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.txt'];

    else
        if regexpbl(ghcnType,'cdf')
            %Open file to write to:
            dirStats = ['GHCN_cdf_stats_' ghcnType];
            dirReport = fullfile(dirTs{kk}, dirStats);
            if ~exist(dirReport,'dir')
                mkdir(dirReport);
            end
            fileGeoStats = ['GHCN_cdf_stn_geoStats_' metVar '_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.asc'];
            nameReport = ['GHCN_cdf_tableStats_' metVar '_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.txt'];

            %Decide on number of bins:
            nBins = 100;
            for jj = 1 : 8
                nBins = nBins - 10;
                suffData = nan( numel(gridReform{kk}(:,1)),1);
                for ii = 1 : numel(gridReform{kk}(:,1))
                    if numel(find(~isnan(gridReform{kk}{ii,5}(:,2:13)) == 1)) > nBins
                        suffData(ii) = 1;
                    end
                end
                if numel(find(~isnan(suffData) == 1)) > 0.75*numel(gridReform{kk}(:,1))
                    nBins = nBins - 10;
                    break
                end
                if jj == 8 && numel(find(~isnan(suffData) == 1)) == 0
                    warning('GHCN_2_grid_cmpr:needMoreData',['There are no '...
                        'stations with sufficient data to do CDF analysis '...
                        'for the current gridded data. Other analysis may still be conducted.']);
                    return
                end
            end

            if nBins > 50
                nBins = 50;
            end
            
            if nBins < 10
               error('GHCN_2_grid_cmpr:nBins',['The number of bins is ' num2str(nBins) ', which should not be possible.']); 
            else
                disp(['The number of bins being used in CDf analysis is ' num2str(nBins) '.']);
            end

            stnStats  =  stn_cdf_stats(gridReform{kk}, stnRec{kk}, nBins, metVar, dirReport);
            mnthStats = mnth_cdf_stats(gridReform{kk}, stnRec{kk}, nBins, metVar, dirReport);
            aggStats  =  agg_cdf_stats(gridReform{kk}, stnRec{kk}, nBins, metVar, dirReport);

            strType = 'CDFs';

            %TAYLOR DIAGRAMS ONLY APPLICABLE FOR COMPARING MULTIPLE GCMS OR METHODS:
            % taylorStats = taylor_stats(gridData, stnRec{kk}, nBins, metVar, dirReport);
            % [hp ht axl] = taylordiag(STDs,RMSs,CORs,['option',value]);

        else
            %Open file to write to:
            dirStats = ['GHCN_ts_stats_' ghcnType];
            dirReport = fullfile(dirTs{kk}, dirStats);
            if ~exist(dirReport,'dir')
                mkdir(dirReport);
            end


            %Calculate ERROR:
            aggStats  =  agg_stats_calc_v2(gridReform{kk}, stnRec{kk});
            mnthStats = mnth_stats_calc_v2(gridReform{kk}, stnRec{kk});
            stnStats  =  stn_stats_calc_v2(gridReform{kk}, stnRec{kk});

            strType = 'time-series';
            fileGeoStats = ['GHCN_ts_stn_geoStats_' metVar '_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.asc'];
            nameReport = ['GHCN_ts_tableStats_' metVar '_' ghcnType '_' num2str(yrs(1)) '_' num2str(yrs(2)) '.txt'];
        end
    end
    
    %WRITE STATION STATS TO CSV FILE FOR USE IN ARCMAP:
    disp([char(10) 'The GHCN to downscaled data comparison script will ' ...
        'write two files with statistics calculated on ' strType ' from ' dataConMessage ...
        ' GHCN data: ' char(10) ...
        'The .txt file stores all computed stats in a format readable in '...
        'a text editor' char(10) 'The .asc file contains statistics binned '...
        'by GHCN station in a delimited format readable by ArcMap.' char(10)]);

    %Use uniform format
    statSpec = '%09.5f';
    %Open file to write to:
    pathGeoStats = fullfile(dirReport, fileGeoStats);

    print_geostats_v2(pathGeoStats, stnRec{kk}, stnYrs{kk}, stnStats);

    pathReport = fullfile(dirReport, nameReport);
    fidReport = fopen(pathReport, 'w');

    %write notes and header to all-stats file
    repIntroStr = ['This report provides Mean Error (bias), ' ...
        'Mean Absolute Error (MAE), Mean Absolute Percent Error (MAPE), ' ...
        'Root Mean Square Error (RMSE), ' ...
        'Weighted Mean Absolute Percent Error (WMAPE), and ' ... 
        'Nash-Sutcliffe Efficiency (NSE) ' ...
        'statistics computed on CDFs from the downscaled time-series data in folder ' dirTs{kk} ...
        ' and the Global Historical Climatological Network (GHCN).' char(10) char(10) ...
        'The range of years included in the statistics is ' num2str(yrs(1)) ...
        ' to ' num2str(yrs(2)) '.' char(10)];
    fprintf(fidReport, '%s', [repIntroStr char(10)]);

    %Metric Definitions and units:
    repFormulaStr = ['Statistical Formulas:' char(10) ...
        'Bias = (1/n)*SUM(M_t - O_t)' char(10) ...
        'MAE = (1/n)*SUM(|M_t - O_t|)' char(10) ...
        'MAPE = (100/n)*SUM(|(M_t - O_t)/O_t|)' char(10) ...
        'WMAPE = (100*O_t)*(SUM|(M_t - O_t )/O_t|/SUM(O_t))' char(10) ...
        'RMSE = sqrt(SUM((M_t - O_t)^2)/n)' char(10) ... 
        'NSE = 1 - SUM((M_t - O_t)^2)/SUM((O_t - O_avg)^2)' char(10) ...
        'where O_t is the GHCN data and F_t is the gridded (downscaled) data.' char(10)];
    fprintf(fidReport, '%s', [repFormulaStr char(10)]);

    repUnitStr = ['Units:' char(10) ...
        '[Bias] = ' units ',' char(10) ...
        '[MAE] = '  units ',' char(10) ...
        '[MAPE] = unitless (percent),' char(10) ...
        '[WMAPE] = unitless (percent),' char(10) ...
        '[RMSE] = ' units ', and' char(10) ...
        '[NSE] = unitless (negative infinity to 1).' char(10)];
    fprintf(fidReport, '%s', [repUnitStr char(10)]);

    %Print GHCN ID of intentionally excluded stations to report:
    print_exclude_ghcn(fidReport, stnEx);

    %Print GHCN ID of stations where gridded data is exclusively NaN to report:
    print_del_ghcn(fidReport, recStnNan,'nan');

    if regexpbl(ghcnType,{'cdf','joint'})
         fprintf(fidReport, '%s', [num2str(nBins) ' bins were used in analysis.' char(10) char(10)]);
         
       	if regexpbl(ghcnType,'joint')
            fprintf(fidReport, '%s', ['Adjusted temperature data and original precipitation GHCN data were used.' char(10) char(10)]);
            
            fprintf(fidReport, '%s', ['Statistics on the correlation between precipitation and mean temperature:' char(10)]);
            %Print Aggregated Pre-Tmp Correlations to file:
  
            agg_stats_print_v2(fidReport, nCorrAggCom, statSpec, length(gridReform{kk}(:,1)) );
            mnth_stats_print_v2(fidReport, nCorrMnthCom, statSpec);
            stn_stats_print_v2(fidReport, gridReform{kk}(:,1), nCorrStnCom, statSpec);
%             agg_stats_print_v2(fidReport, nCorrAggGHCN, statSpec, length(gridReform{kk}(:,1)) );
%             agg_stats_print_v2(fidReport, nCorrAggGrid, statSpec, length(gridReform{kk}(:,1)) );
%             %Print Monthly pre-Tmp Correlations to file:
%             mnth_stats_print_v2(fidReport, nCorrMnthGHCN, statSpec);
%             mnth_stats_print_v2(fidReport, nCorrMnthGrid, statSpec);
%             %Print Station Pre-Tmp Correlations to file:
%             stn_stats_print_v2(fidReport, gridReform{kk}(:,1), nCorrStnGHCN, statSpec);
%             stn_stats_print_v2(fidReport, gridReform{kk}(:,1), nCorrStnGrid, statSpec);
            
            fprintf(fidReport, '%s', [char(10) char(10)]);
            fprintf(fidReport, '%s', ['Statistics on the differences between GHCN and gridded histograms using ' num2str(nBins) ' bins:' char(10)]);
%             sprintf(fidReport,['The correlation between precipitation and temperature in the aggregated (over time-series and stations) GHCN dataset is ' num2str(nCorrAggGHCN) '.' char(10)]);
%             sprintf(fidReport,['The correlation between precipitation and temperature in the aggregated (over time-series and stations) gridded dataset is ' num2str(nCorrAggGrid) '.' char(10)]);
        end
    end
    
    %Print Aggregated Error to file:
    agg_stats_print_v2(fidReport, aggStats, statSpec, length(gridReform{kk}(:,1)) );
    %Print Monthly Statistics to file:
    mnth_stats_print_v2(fidReport, mnthStats, statSpec);
    %Print Station Stats to file:
    stn_stats_print_v2(fidReport, gridReform{kk}(:,1), stnStats, statSpec);

    fclose(fidReport);
end

