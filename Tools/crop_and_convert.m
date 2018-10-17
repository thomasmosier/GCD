pathIn = '/Users/thomas/Library/Mobile Documents/com~apple~CloudDocs/climate_data_temp/ERA_daily/rsds/rsds_day_ERAi_reanalysis_19790101-20171231.nc';
varLd = 'rsds';
yrsLd = [1986,2015];

lonLd = [63, 71];
latLd = [34, 39];

unitConv = 10^(-6); %Convert W per m^2 to MJ
unitsOut = 'MJ m**-2';

prec = 1;
wrtTyp = 'netcdf';


%READ DATA
sData = read_geodata_v2(pathIn, varLd, lonLd, latLd, yrsLd, 1, 'in');


%CONVERT UNITS
varAtt = ['att' varLd];
varUnits = 'Units';
sData.(varLd) = sData.(varLd)*unitConv;
unitsIn = find_att(sData.(varAtt), varUnits);
disp(['The input units are ' unitsIn '. These will be set to ' unitsOut ...
    ' through a multiplicative factor of ' num2str(unitConv) '.']);
sData.(varAtt) = set_att(sData.(varAtt), varUnits, unitsOut);


%WRITE OUTPUT
[foldIn, fileIn, extIn] = fileparts(pathIn);
foldOut = fullfile(foldIn, 'new_units');
if ~exist(foldOut, 'dir')
   mkdir(foldOut); 
end
pathOut = fullfile(foldOut, [fileIn, extIn]);
if exist(pathOut, 'file')
   delete(pathOut); 
end
write_geodata_v2(pathOut, sData, prec, wrtTyp, 'var', varLd);

