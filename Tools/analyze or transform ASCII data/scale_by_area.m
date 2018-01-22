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


scl = 1000; %Multiplies the data by this factor when processing.
fun = 'd'; %Multiply input grid by area ('m') or divide grid by area ('d');
typeGet = 'file'; %Select individual 'file' or 'folder' of data.

strSlct = ['Select the ESRI ASCII ' typeGet ' to scale by the area of each cell.'];
if ~isempty(regexpi(typeGet,'file'))
    [fileTemp, path] = uigetfile({'*.asc';'*.txt'}, strSlct);
    file = cell(1,1);
    file{1} = fileTemp;
elseif ~isempty(regexpi(typeGet,'fold'))
    fold = uigetdir(pwd,strSlct);
    
    fileTemp = dir(fullfile(fold,'*.asc'));
    file = struct2cell(fileTemp);
    file = file(1,:);
else
    error('scale_by_area:typeGet',[char(39) typeGet char(39) ' is not a recognized option.']);
end

for ii = 1 : numel(file)
    file{ii} = fullfile(path,file{ii});
    disp([char(39) file{ii} char(39) ' has been selected as the file to scale by its area.']);

    [data, hdr, meta] = read_ESRI(file{ii});

    [lat, lon] = ESRI_hdr2geo(hdr, meta);
    area = area_geodata(lon, lat);

    if ~isempty(regexpi(fun,'m'))
        data = data .* area;
    elseif  ~isempty(regexpi(fun,'d'))
        data = data ./ area;
    else
        error('scale_area:unknown',ln80(['The string option ' char(39) fun ...
            char(39) ' has not been programmed for.']));
    end

    ind = regexpi(file{ii},'\.');
    fileNew = [file{ii}(1:ind(end)-1) '_' fun 'area.asc']; 

    disp(ln80(['The area-scaled data are being written to ' char(39) ...
        fileNew char(39) '.']));
    if scl ~= 1
        disp(['The data have been multipled by ' num2str(scl) '.']);
        data = data * scl;
    end

    %Identify how many decimals to write:
    dec = dec_write(data);
    disp([num2str(dec) ' decimal places are being recorded in the output data.']);

    write_ESRI_v3(data, hdr, fileNew, dec);
end