function dataModBcOut = bc1d_day(dataHistRef, datesHistRef, dataHistMod, datesHistMod, dataMod2Bc, datesMod2Bc, nDays, type)

%Floor all dates to ensure they are same type and comparable
datesHistRef = floor(datesHistRef);
datesHistMod = floor(datesHistMod);
datesMod2Bc = floor(datesMod2Bc);

blWarn = 1;
[~, msgid] = lastwarn;
if regexpbl(msgid, 'bc1dDay')
    blWarn = 0;
end

%Check that data are daily:
if (numel(datesHistRef(1,:)) < 3 || nanmean(diff(datesHistRef(:,3))) > 4) &&  blWarn == 1
    error('bc1dDay:timeResHistRef', ['The time resolution of the historic ' ...
        'reference time-series appears to not be daily. Daily data are needed ' ...
        'for this bias correction method.']);
end
if (numel(datesHistMod(1,:)) < 3 || nanmean(diff(datesHistMod(:,3))) > 4) &&  blWarn == 1
    error('bc1dDay:timeResHistMod', ['The time resolution of the historic ' ...
        'simulation time-series appears to not be daily. Daily data are needed ' ...
        'for this bias correction method.']);
end
if (numel(datesMod2Bc(1,:)) < 3 || nanmean(diff(datesMod2Bc(:,3))) > 4) &&  blWarn == 1
    error('bc1dDay:timeResHistMod', ['The time resolution of the ' ...
        'simulation time-series to bias correct appears to not be daily. ' ...
        ' Daily data are needed for this bias correction method.']);
end


%Define unique quantile mapping for each 31 day window (Loop over each day of the year)
%Find unique days of the year (fur creating QM transfer function)
datesMod2BcLp = floor(unique(datesMod2Bc(:,2:3), 'rows'));

%Initialize bias corrected data:
dataModBcOut = nan([numel(datesMod2Bc(:,1)), 1]);

%Loop over each unique month-day to create quantile mapping
for nn = 1 : numel(datesMod2BcLp(:,1))
    %Find current 31 day window (index refers to input dates in for-loop):
    indTCurr = nn + (-nDays : nDays);

    %Wrap edges of 31 day window:
    indTNeg = find(indTCurr < 1);
    indTPos = find(indTCurr > numel(datesMod2BcLp(:,1)));

    indTCurr(indTNeg) = numel(datesMod2BcLp(:,1)) + indTCurr(indTNeg);
    indTCurr(indTPos) = indTCurr(indTPos) - numel(datesMod2BcLp(:,1));

    %Find current QM dates:
    dates2BcCurr = datesMod2BcLp(indTCurr, :);

    %Collect time-series indices corresponding to QM dates:
    indHistModCurr = find(ismember(datesHistMod(:,2:3),        dates2BcCurr, 'rows') == 1);
    indHistRefCurr = find(ismember(datesHistRef(:,2:3),        dates2BcCurr, 'rows') == 1);
    indMod2BcCurr  = find(ismember( datesMod2Bc(:,2:3), datesMod2BcLp(nn,:), 'rows') == 1);

    if (numel(indHistModCurr) < 30*28 && nn == 1) &&  blWarn == 1
        warning('bc1dDay:lowModTs', ...
            [num2str(numel(indHistModCurr)) ' time-series elements ' ...
            'were found in the historical reference dataset to use in training the bias correction for the ' ...
            'current time step. Usually 30 years of data is recommended.']);
    end
    if (numel(indHistRefCurr) < 30*28 && nn == 1) &&  blWarn == 1
        warning('bc1dDay:lowRefTs', ...
            [num2str(numel(indHistRefCurr)) ' time-series elements ' ...
            'were found in the historical GCM simulation to use in training the bias correction for the ' ...
            'current time step. Usually 30 years of data is recommended.']);
    end

    if any(~isnan(dataModBcOut(indMod2BcCurr))) &&  blWarn == 1
        error('bc1dDay:notNanOutTs',['Non-nan values are '...
        'present in GCM projection simulation before bias correction. This is unexpected.']);
    else
        if regexpbl(type, 'eQM')
            dataModBcOut(indMod2BcCurr) = e_qm(dataHistRef(indHistRefCurr), dataHistMod(indHistModCurr), ...
                dataMod2Bc(indMod2BcCurr), type);
        elseif regexpbl(type, 'q2q')
            [dataModBcOut(indMod2BcCurr)] = e_q2q(dataHistRef(indHistRefCurr), dataHistMod(indHistModCurr), ...
                dataMod2Bc(indMod2BcCurr));
        end
    end
end
clear nn

if (any2d(isnan(dataModBcOut)) && ~all(isnan(dataModBcOut))) &&  blWarn == 1
    error('bc1dDay:nanOutQm',['Nan values are '...
        'present in GCM projection simulation bias corrected with QM. This is unexpected.']);
end