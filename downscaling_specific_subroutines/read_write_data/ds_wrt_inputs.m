function ds_wrt_inputs(sTVar, sDs, sPath)     


if numel(sDs.fldsTVar(:)) ~= numel(sTVar(:))
    error('dsWrtInputs:diffNumberInputs', [num2str(numel(sTVar(:))) ...
        ' data structures were input but only ' num2str(numel(sDs.fldsTVar(:))) ...
        ' field descriptions.']);
end

for kk = 1 : numel(sDs.fldsTVar(:))
    if regexpbl(sDs.fldsTVar{kk}, sDs.wrtOut)
        foldCurr = fullfile(sPath.output, 'extras', sDs.fldsTVar{kk});
        if ~exist(foldCurr, 'dir')
            mkdir(foldCurr);
        end
        
        nmOut = ds_file_nm(sPath, sDs, kk);
        
        fileCurr = [sDs.varDs, '_', nmOut, '_', ...
            sDs.region];

        if isfield(sTVar{kk}, sDs.varDs)
            write_geodata_v2(fullfile(foldCurr, fileCurr), sTVar{kk}, sDs.wrtPrec, sDs.wrtFmt, 'var', sDs.varDs);
        end
    end
end
clear kk