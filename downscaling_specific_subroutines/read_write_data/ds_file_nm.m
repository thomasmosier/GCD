function nmOut = ds_file_nm(sPath, sDs, indFld)

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