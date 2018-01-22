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

function sSimPBc = bc1d_geodata(sSim2Bc, sSimH, sRefH, var, yrs, type)

varLat = 'latitude';
varLon = 'longitude';
varDate = 'date';

%Create empirical cumulative distribution functions (ECDFs) at same spatial 
%resolution for 'sGcm' and 'sObs'.  Both datasets are necessary as inputs
%to ensure they are at same spatial resolution.

% %Identify which grid input is higher resolution.  Resample to the lower
% %resolution.
% if (~isequal(sSimH.(varLat), sRefH.(varLat)) || ~isequal(sSimH.(varLon), sRefH.(varLon))) && all([numel(sSimH.(varLon)), numel(sSimH.(varLat)), numel(sRefH.(varLon)), numel(sRefH.(varLat))] ~= 1)
%     warning('q_map_geodata:resample',['The input grids have different '...
%         'resolutions so the higher resolution grid is being resampled to the lower resolution.']);
%     [sSimH, sRefH, nRe] = upscale_geodata_v2(sSimH, sRefH, var); %Last argument states that if there is a nan in a sum of numbers during re-sampling, the nan wont affect the output
%     if nRe == 1
%         sSimH.(var) = sSimH.(var);
%         sSimH.(varLat) = sSimH.([varLat 'Re']);
%         sSimH.(varLon) = sSimH.([varLon 'Re']);
%     elseif nRe == 2
%         sRefH.(var) = sRefH.(var); 
%         sRefH.(varLat) = sRefH.([varLat 'Re']);
%         sRefH.(varLon) = sRefH.([varLon 'Re']);
%     end
% end

%Find lengths of data:
nLat = length(sSim2Bc.(varLat));
nLon = length(sSim2Bc.(varLon));


%Initialize bias-corrected data structure array:
sSimPBc = sSim2Bc;
if isfield(sSimPBc,'val')
    sSimPBc = rmfield(sSimPBc,'val');
end
if isfield(sSimPBc,'cdf')
    sSimPBc = rmfield(sSimPBc,'cdf');
end
sSimPBc.(var) = nan(size(sSim2Bc.(var)),'single');


%Find time-series elements to compare during bias correction:
%For historic model:
indSimH = find(sSimH.(varDate)(:,1) >= min(yrs) & sSimH.(varDate)(:,1) <= max(yrs));
indRefH = find(sRefH.(varDate)(:,1) >= min(yrs) & sRefH.(varDate)(:,1) <= max(yrs));

if numel(indRefH) ~= numel(indSimH)
    [indSameSimH, ~] = ismember(sSimH.(varDate)(:,1), sRefH.(varDate)(:,1), 'rows');
    [indSameRefH, ~] = ismember(sRefH.(varDate)(:,1), sSimH.(varDate)(:,1), 'rows');
    
    indSimH = intersect(indSimH, sort(find(indSameSimH ~= 0)));
    indRefH = intersect(indRefH, sort(find(indSameRefH ~= 0)));
	warning('eQMgeodata:unequalYrs', ['The number of years present in the '...
       'two datasets being compared are not equal. They are going to be set to the common years of ' num2str(min(sSimH.(varDate)(indSimH,1))) '-' num2str(max(sSimH.(varDate)(indSimH,1)))]); 
end


%Do quantile mapping bias-correction
for ii = 1 : nLat
	for jj = 1 : nLon
        %Find indices for each output grid cell:
        [~, indLatSimH] = min(abs(sSim2Bc.(varLat)(ii) - sSimH.(varLat)));
        [~, indLonSimH] = min(abs(sSim2Bc.(varLon)(jj) - sSimH.(varLon)));
        
        [~, indLatRefH] = min(abs(sSim2Bc.(varLat)(ii) - sRefH.(varLat)));
        [~, indLonRefH] = min(abs(sSim2Bc.(varLon)(jj) - sRefH.(varLon)));
        
        if regexpbl(type, 'month')
            sSimPBc.(var)(:,ii,jj) = bc1d_month(...
                sRefH.(var)(indRefH,indLatRefH,indLonRefH), sRefH.(varDate)(indRefH,:), ...
                sSimH.(var)(indSimH,indLatSimH,indLonSimH), sSimH.(varDate)(indSimH,:), ...
                sSim2Bc.(var)(:,ii,jj), sSim2Bc.(varDate),...
                type);
        elseif regexpbl(type, {'day','window'})
            %Find number of days (i.e. +/- current)
            nDyWin = regexpi(type,'[\d]{1,3}day','match');
%             nDyWin = regexp(type,'-?\d+\.?\d*|-?\d*\.?\d+','match');
            if isempty(nDyWin) 
                nDyWin = 15;
            elseif iscell(nDyWin)
                if numel(nDyWin) > 1
                    error('bc1d:multDayWindows',[num2str(numel(nDyWin)) ' day window indicators were found.'])
                else
                    nDyWin = str2double(char(regexpi(nDyWin{1},'[\d]{1,3}','match')));
                end
            end

            sSimPBc.(var)(:,ii,jj) = bc1d_day(...
                sRefH.(var)(indRefH,indLatRefH,indLonRefH), sRefH.(varDate)(indRefH,:), ...
                sSimH.(var)(indSimH,indLatSimH,indLonSimH), sSimH.(varDate)(indSimH,:), ...
                sSim2Bc.(var)(:,ii,jj), sSim2Bc.(varDate),...
                nDyWin, type);
        else
            error('eQmGeodata:unknownType',['eQM type ' type ' has not been programmed for.']);
        end
	end
end

%Various plots for testing:
% plot(linspace(0,1,nMnCDF),sObs.val{3,5},'s',linspace(0,1,nMnCDF),sGcm.val{3,5},'s')
% plot(linspace(0,1,nMnCDF),map{3,5},'s')
% plot(linspace(0,1,nMnCDF),map{3,5},'s')