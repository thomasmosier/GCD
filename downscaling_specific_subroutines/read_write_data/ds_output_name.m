function nmOut = ds_output_name(sDs, strData, mnth, typOut)




if mnth < 10
    strMnth = ['0' num2str(mnth)];
else
    strMnth = num2str(mnth);
end
    

%Create date string:
if strcmpi(typOut, 'ds')
    if regexpbl(sDs.wrtFmt, 'nc')
        fileDate = {[ '_' ...
            num2str(min(sDs.yrsDs)) strMnth '01' '-' ...
            num2str(max(sDs.yrsDs)) strMnth num2str(eomday(max(sDs.yrsDs), mnth))]};
%     elseif regexpbl(sDs.wrtFmt, 'asc')
%         yrs = (min(sDs.yrsDs):max(sDs.yrsDs));
%         fileDate = cell(numel(yrs), 1);
% 
%         for ii = 1 : numel(fileDate)
%             fileDate{ii} = [num2str(yrs(ii)) '_' strMnth];
%         end
    else
        fileDate = {''};
    end
    
    tStep = sDs.timestep;
elseif strcmpi(typOut, 'clim')
    tStep = 'clim';
    
    fileDate = {[num2str(min(sDs.yrsDs)) 'thru' num2str(max(sDs.yrsDs)) '_' strMnth]};
else
    error('dsOutputName:typeUnknown',['The output type ' typOut ' has not been programmed for in this function.']);
end

nmRt = [sDs.varDs '_' tStep '_' strData '-ds-' sDs.region '_' sDs.method '_intrp-' sDs.intrp];

nmOut = cell(numel(fileDate), 1);
for ii = 1 : numel(fileDate)
    if isempty(fileDate{ii})
        nmOut{ii} = nmRt;
    else  
        nmOut{ii} = [nmRt, '_' fileDate{ii}];
    end
end