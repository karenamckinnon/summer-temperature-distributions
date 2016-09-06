function [ISD QR] = mainHourly(yearStart, yearEnd, months, percentiles, varname, stationList, ...
	timesToSample, timeForQR, runNumber, nIntervals, njitter, nboot)

% Sample hourly data at different measurement times and determine effect on QR slopes,
% as compared to the warming / variability signal
 
% user needs to fill these in 
csvDir = ''; % location of csv files with ISD data
cacheDir = '';
figDir = '';

fid = fopen(stationList);
headers = textscan(fid, '%q', 11, 'delimiter', ',');
vals = textscan(fid, '%d %d %s %s %s %s %f %f %f %d %d', 'delimiter', ',');
fclose(fid);
allUSAF = reshape(sprintf('%06d', vals{1}),[6 length(vals{1})])';
allWBAN = reshape(sprintf('%05d', vals{2}),[5 length(vals{2})])';
 
startInd = 1;

opts.yearStart = yearStart;
opts.yearEnd = yearEnd;
opts.csvDir = csvDir;
opts.cacheDir = cacheDir;
opts.startInd = startInd;
opts.months = months;

% sample hourly data at two different times
% (will take a long time if not previously calculated)

[ISD] = sampleHourlyWrapper(allUSAF, allWBAN, timesToSample, runNumber, nIntervals, opts);

ISD0 = ISD;

% subselect one obs time for the rest of the analysis
ISD.TMAX = squeeze(ISD0.TMAX(:, timesToSample == timeForQR, :));
ISD.TMIN = squeeze(ISD0.TMIN(:, timesToSample == timeForQR, :));

ISD = doQC(ISD, months, varname);

% Do QR
QR = QRwrapper(ISD, varname, nIntervals, months, runNumber, percentiles, ...
	njitter, cacheDir, nboot);

return

