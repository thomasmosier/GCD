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

function [corrMat] = shep_weight_array(stnData, gridData, gridHdr, date, p, r)


%This function uses Shepard's weighting (see 1968 paper) to produce a
%correction matrix at each time-step, where the list of correction factors
%is the difference between GHCN station records and downscaled data
%surface.

%Changes that need to be made:
mnth = date(1);
yr = date(2);
    %Integrate with ''
    
[lonGrid, latGrid] = geo_mesh(gridHdr);

%stnData and gridData are cell arrays with the form:
    %{i,1} is the GHCN station ID,
    %{i,2} is the cell-center latitude, 
    %{i,3} is the cell-center longitude,
    %{i,4} is the station elevation, and 
    %{i,5} is an m by 13 cell-array with yearly data (first column is year).

% %Determine first and last years present in station records:
% for ii = 1 : length( gridData(:,5) )
%     if ii == 1
%         startYr = gridData{ii,5}(1,1);
%             startInd = ii;
%         endYr = gridData{ii,5}(end,1);
%             endInd = ii;
%     else
%         if gridData{ii,5}(1,1) < startYr
%             startYr = gridData{ii,5}(1,1);
%                 startInd = ii;
%         end
%         if gridData{ii,5}(end,1) > endYr
%             endYr = gridData{ii,5}(end,1);
%                 endInd = ii;
%         end
%     end
% end

% yrs = startYr : 1 : endYr;

% %Initialize corrAry:
% corrAry = cell( length(yrs), length(mnths) );
% [ corrAry{:,:} ] = deal( zeros(length(gridAry(:,1)),3) );
% corrMeta = zeros(length(yrs)*length(mnths), 3);

% %Loop over years:
% for ii = 1 : length(yrs)
% 
%     %Loop over months:
%     for jj = 1 : length(mnths)
        
%         %Enter date into metadata matrix:
%         corrMeta((ii-1)*12+jj, 1) =   yrs(ii);
%         corrMeta((ii-1)*12+jj, 2) = mnths(jj);
        
        %Find stations with meteorological values for the current year:
        
        %Initialize:
        iiGridMat = nan(length(gridData(:,1)),3);
        iiStnMat  = nan(length(gridData(:,1)),3);
        
        useCntr = 1;    %Cntr for loop over meteorological stations
        
        for kk = 1 : length(gridData(:,1))
            currYr = find( gridData{kk,1}(:,1) == yr );
            if ~isempty(currYr) && ~isnan(gridData{kk,1}(currYr, mnth+1))
                %Copy latitudes:
                iiGridMat(useCntr,1) = gridData(kk,2);
                 iiStnMat(useCntr,1) =  stnData(kk,2);

                %Copy longitudes:
                iiGridMat(useCntr,2) = gridData(kk,3);
                 iiStnMat(useCntr,2) =  stnData(kk,3);  
                 
                %Copy meteorological values:
                iiGridMat(useCntr,3) = gridData{kk,5}(currYr, mnth+1);
                 iiStnMat(useCntr,3) =  stnData{kk,5}(currYr, mnth+1);
                 
                useCntr = useCntr + 1;
            end
        end
        
        %delete empty rows:
        if useCntr <= length(iiStnMat(:,1))
            iiGridMat(useCntr,:) = [];
             iiStnMat(useCntr,:) = [];
        end
   
        %Initialize correction matrix
        corrMat = nan( gridHdr(2),gridHdr(1) );
    
        if ~isempty( iiStnMat(:,1) )    %Only proceed if at least one station exists for current timestep.
            for kk = 1 : length( corrMat(:,1) ) %interate over rows of corrMat
                for ll = 1 : length( corrMat(1,:) ) %iterate over columns of corrMat

                    %Check if there is a station within the current grid cell 
                    stnInCell = intersect( ...
                        find( iiGridMat(:,1)==latGrid(kk,1) ), ...
                        find( iiGridMat(:,2)==lonGrid(1,ll) ) ...
                        );

                    if ~isempty(stnInCell)  %If grid pt contains stn, correction factor is station bias: 
                        corrMat(kk,ll) = iiStnMat(stnInCell, 3) - iiGridMat(stnInCell, 3);
                        
                    else %If grid pt does not contain a station, compute Shepard's weight:
                        deltaLat = iiStnMat(:,1) - latGrid(kk,1);   %Vector
                        deltaLon = iiStnMat(:,2) - lonGrid(1,ll);   %Vector
                        
                        %Vector of distances between current cell and each
                        %station with value at current time step.
                        dJK = ...
                            r*2*asin(sqrt( sind(deltaLat/2).^2 ...
                            + cosd(latGrid(kk,1)).*cosd(iiStnMat(:,1)).*sind(deltaLon/2).^2 ));

                        corrMat(kk,ll) ...
                            = sum(dJK.^(-p)*(iiStnMat(:,3) - iiGridMat(:,3)))...
                            /sum(dJK.^(-p));
                    end

                end
            end

%             %Meta matrix gets '1' if a correction has been applied to the
%             %time step.
%             corrMeta((ii-1)*12+jj, 3) = 1;
        else    %Case where no stn data exists for current time step.
%             %Meta matrix gets '0'
%             corrMeta((ii-1)*12+jj, 3) = 0;
        end
        
%         %Write correction Matrix for current time-step to corrAry
%         corrAry{ii,jj} = corrMat;
    

    
        
%     end     %End of month loop
% end     %End of year loop
   


end     %End of function!