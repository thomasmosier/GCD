function nmOut = ds_file_nm(sPath, sDs, indFld, varargin)

date = [];
if ~isempty(varargin(:))
    for ii = 1 : numel(varargin(:))
        if strcmpi(varargin{ii}, 'date')
            date = varargin{ii+1};
        end
    end
end

if isfield(sPath, ['nm_' sDs.fldsTVar{indFld}])
    if iscell(sPath.(['nm_' sDs.fldsTVar{indFld}]))
        if numel(sPath.(['nm_' sDs.fldsTVar{indFld}])(:,1)) == 1
            nmOut = sPath.(['nm_' sDs.fldsTVar{indFld}]){1,1};
        else
            nmOut = sPath.(['nm_' sDs.fldsTVar{indFld}]){sDs.indSim}{1};
        end
    elseif ischar(sPath.(['nm_' sDs.fldsTVar{indFld}]))
        nmOut = sPath.(['nm_' sDs.fldsTVar{indFld}]);
    end
else
    indUnd = regexpi(sDs.fldsTVar{indFld}, '_');
    if regexpbl(sDs.fldsTVar{indFld}(indUnd(end):end), 'bc')
        testFld = ['nm_' sDs.fldsTVar{indFld}(1:indUnd(end)-1)];
    else
        testFld = ['nm_' sDs.fldsTVar{indFld}];
    end
    
    if isfield(sPath, testFld);
        if iscell(sPath.(testFld))
            if numel(sPath.(testFld)(:,1)) == 1
                nmOut = sPath.(testFld){1,1};
            else
                nmOut = sPath.(testFld){sDs.indSim}{1};
            end
        elseif ischar(sPath.(testFld))
            nmOut = sPath.(testFld);
        end
        nmOut = [nmOut sDs.fldsTVar{indFld}(indUnd(end):end)];
    else
        nmOut = sDs.fldsTVar{indFld};
    end
end

if ~isempty(date)
    date = floor(date);
    indDate4 = regexpi(nmOut, '\d{4}');
    indDate8 = regexpi(nmOut, '\d{8}');
    if ~isempty(indDate8) && numel(indDate8) == 2
        yr1 = num2str(min(date(:,1)));
        mnth1 = num2str(min(date(:,2)));
        if numel(mnth1) == 1
            mnth1 = ['0' mnth1];
        end
        day1 = num2str(min(date(:,3)));
        if numel(day1) == 1
            day1 = ['0' day1];
        end
        date1 = [yr1 mnth1 day1];
        
        yr2 = num2str(max(date(:,1)));
        mnth2 = num2str(max(date(:,2)));
        if numel(mnth2) == 1
            mnth2 = ['0' mnth2];
        end
        day2 = num2str(max(date(:,3)));
        if numel(day2) == 1
            day2 = ['0' day2];
        end
        date2 = [yr2 mnth2 day2];
        
        nmOut(indDate8(1):indDate8(1)+7) = date1;
        nmOut(indDate8(2):indDate8(2)+7) = date2;
    elseif isempty(indDate8) && ~isempty(indDate4)
        warning('dsFileNm:monthlyTimestep','The data are likely in monthly timestep naming format. This clause has been reserved but not filled in.')
    end
end