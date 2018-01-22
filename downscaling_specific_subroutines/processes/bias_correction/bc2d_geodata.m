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

function  [sSimPBc, sSimPBcV] = bc2d_geodata(sSim2Bc, sSim2BcV, sSimH, sSimHV, sRefH, sRefHV, var, varV, yrsBase, type)


%Create empirical cumulative distribution functions (ECDFs) at same spatial 
%resolution for 'sGcm' and 'sObs'.  Both datasets are necessary as inputs
%to ensure they are at same spatial resolution.

%Inputs:
%'sModV1' = first variable of the structure array (Downscaling Package 
    %format) that map is created from
%'sObsV1' = first variable of the structure array (Downscaling Package 
    %format) that map is created to
%'sModV2' = second variable of the structure array (Downscaling Package 
    %format) that map is created from
%'sObsV2' = second variable of the structure array (Downscaling Package 
    %format) that map is created to
%'nBins' = number of bins to use in producing empirical cumulative 
    %distribution function 
%'window' = Indicates the size of the window used when creating the map 
    %(Either '1' or '3')
%'varargin' = if 1, scales s2 to s1 grid; if 2, scales s1 to s2 grid;
    %if empty, scales to higher-spatial resolution grid

varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';

%Identify which grid input is higher resolution. Resample to the lower
%resolution.
flagRe = 0;
if ~isequal(sSimHV.(varLat), sRefHV.(varLat)) || ~isequal(sSimHV.(varLon), sRefHV.(varLon))
    [sSimHV, sRefHV, nRe] = upscale_geodata_v2(sSimHV, sRefHV, varV); %Last argument states that if there is a nan in a sum of numbers during re-sampling, the nan wont affect the output
    flagRe = 1;
end
if ~isequal(sSimH.(varLat), sRefH.(varLat)) || ~isequal(sSimH.(varLon), sRefH.(varLon))
    [sSimH, sRefH, ~] = upscale_geodata_v2(sSimH, sRefH, var, nRe);
    flagRe = 1;
end
if flagRe == 1
    warning('eJBCGeodata:resample',['The input grids have different '...
        'resolutions so the higher resolution grid is being resampled to the lower resolution.']);
    if nRe == 1
        if isfield(sSimHV, [varV 'Re'])
            sSimHV.(varV) = sSimHV.([varV 'Re']);
            sSimHV.(varLat)  = sSimHV.([varLat 'Re']);
            sSimHV.(varLon)  = sSimHV.([varLon 'Re']);
        end
        
        if isfield(sSimH, [var 'Re'])
            sSimH.(var) = sSimH.([var 'Re']);
            sSimH.(varLat)  = sSimH.([varLat 'Re']);
            sSimH.(varLon)  = sSimH.([varLon 'Re']);
        end
    elseif nRe == 2
        if isfield(sRefHV, [varV 'Re'])
            sRefHV.(varV) = sRefHV.([varV 'Re']);
            sRefHV.(varLat)  = sRefHV.([varLat 'Re']);
            sRefHV.(varLon)  = sRefHV.([varLon 'Re']);
        end
        
        if isfield(sRefH, [var 'Re'])
            sRefH.(var) = sRefH.([var 'Re']);
            sRefH.(varLat)  = sRefH.([varLat 'Re']);
            sRefH.(varLon)  = sRefH.([varLon 'Re']);
        end
    end
end

%Initialize output:
sSimPBc = sSim2Bc;
sSimPBc.(var) = nan(size(sSimPBc.(var)), 'single');

sSimPBcV = sSim2Bc;
sSimPBcV.(var) = nan(size(sSimPBcV.(var)), 'single');


%Find time-series elements to compare during bias correction:
%For historic model:
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

if ~isequal(sSimH.(varDate), sSimHV.(varDate))
    error('eJbcGeodata:diffDateSim',['The historic simulation data for ' ...
        var ' and ' varV ' have different dates. This has not been programmed for.']);
end

if ~isequal(sRefH.(varDate), sRefHV.(varDate))
    error('eJbcGeodata:diffDateRef',['The historic reference data for ' ...
        var ' and ' varV ' have different dates. This has not been programmed for.']);
end

if ~isequal(sSim2Bc.(varDate), sSim2BcV.(varDate))
    error('eJbcGeodata:diffDateRef', ['The projection simulation data for ' ...
        var ' and ' varV ' have different dates. This has not been programmed for.']);
end

nLat = length(sSim2Bc.(varLat));
nLon = length(sSim2Bc.(varLon));

%Construct map (indice 1 = bin from reference variable, indice 2 = correction for all VCor values within the current bin):
for ii = 1 : nLat %Loop over lat points
    for jj = 1 : nLon %Loop over lon points
        %Find indices for each output grid cell:
        [~, indLatSimH] = min(abs(sSim2Bc.(varLat)(ii) - sSimH.(varLat)));
        [~, indLonSimH] = min(abs(sSim2Bc.(varLon)(jj) - sSimH.(varLon)));
        
        [~, indLatRefH] = min(abs(sSim2Bc.(varLat)(ii) - sRefH.(varLat)));
        [~, indLonRefH] = min(abs(sSim2Bc.(varLon)(jj) - sRefH.(varLon)));
        
        if regexpbl(type, 'month')
            mnths = unique(sSim2Bc.(varDate)(:,2));
            for kk = 1 : numel(mnths)
                indRefHCurr = intersect(find(sRefH.(varDate)(:,2) == mnths(kk)), indRefH);
                indSimHCurr = intersect(find(sSimH.(varDate)(:,2) == mnths(kk)), indSimH);
                indSimPCurr = find(sSim2Bc.(varDate)(:,2) == mnths(kk));

                [sSimPBc.(var)(indSimPCurr,ii,jj), sSimPBcV.(varV)(indSimPCurr,ii,jj)] ...
                    = e_JBC(...
                     sRefH.(var)(indRefHCurr,indLatRefH,indLonRefH), sRefHV.(varV)(indRefHCurr,indLatRefH,indLonRefH), ...  
                     sSimH.(var)(indSimHCurr,indLatSimH,indLonSimH), sSimHV.(varV)(indSimHCurr,indLatSimH,indLonSimH), ... 
                     sSim2Bc.(var)(indSimPCurr,ii,jj), sSim2BcV.(varV)(indSimPCurr,ii,jj), type, type);
            end
            clear kk
%         elseif regexpbl(type, {'day','window'})
%               disp('Add functionality...')
%             sSimPBc.(var)(:,ii,jj) = eJBC_day_window(...
%                 sRefH.(var)(indRefH,ii,jj), sRefH.(varDate)(indRefH,:), ...
%                 sSimH.(var)(indSimH,ii,jj), sSimH.(varDate)(indSimH,:), ...
%                 sSimP.(var)(:,ii,jj), sSimP.(varDate),...
%                 type);
        else
            error('eJbcGeodata:unknownType',['eJBC type ' type ' has not been programmed for.']);
        end    
    end
    clear jj
end
clear ii
