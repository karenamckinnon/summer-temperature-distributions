function [GHCND QR] = mainGHCND(yearStart, yearEnd, months, percentiles, varname, ...
	ghcnd_matpath, latRange, lonRange, ... 
	runNumber, njitter, nboot, nIntervals);
% Get GHCND data, and calculate quantile regression slopes
 
% user needs to fill these in
dataDir = ''
cacheDir = '';
figDir = '';  

% get data
GHCND = getGHCND(varname, ghcnd_matpath, yearStart, yearEnd, months, latRange, lonRange, cacheDir);

% Do QC
GHCND = doQC(GHCND, months, varname);

% switch units from 1/10 C to C
GHCND.(varname) = GHCND.(varname)/10;

% Do QR
QR = QRwrapper(GHCND, varname, nIntervals, months, runNumber, percentiles, ...
	njitter, cacheDir, nboot);






