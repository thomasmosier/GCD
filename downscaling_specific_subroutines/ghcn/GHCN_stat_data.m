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

function GHCN_stat_data(category, stats, statSpec, fidStats) 



%This function writes the statistical data in a matrix format to a text
%file.  The rows correpond to the number of items in the category 
%(e.g. n months in a year) and the number of columns equals the number of
%statistics used plus a header column.

%'category' is a vector (potentially cell array, vector of numbers, etc) which label rows.
%'stnStats' is a structure array of the statistical values and statistic
%names.
%'statSpec' defines the format for writing the data (i.e. precision and type)
%'fidStats' is the file ID for the file being written.


nStats = length(stats);

%Loop over the number of data points (rows) for each statistic.
for ii = length(stats(1).data(:)) : -1 : 1

%     if ii == length(stnStats(1).data(:))
%         tempName = 'Avg';
%     else
        tempName = '';
        if iscell(category)
            if ischar(category{ii})
                tempName = category{ii};
            elseif strcmpi(class(category{ii}),'double')
                tempName = num2str(category{ii});
            end
        elseif ischar(category)
            tempName = category;
        elseif strcmpi(class(category),'double')
            tempName = num2str(category(ii));
        else
            warning('GHCN_stat_data:class',['The class of the variable stnNmbr is ' ...
                char(class(category)), ' which does not have a defined conversion option.']);
        end
%     end
    
    strData = [ tempName ' -> '];
    
    %Loop over the number of statistics.
    for jj = 1 : nStats
        strData = [char(strData), num2str(stats(jj).data(ii), statSpec)];
        
        if jj ~= nStats
            strData = [char(strData), ', '];
        end
    end

    fprintf(fidStats, '%s', [char(10) strData]);
end

fprintf(fidStats, '%s', char(10));