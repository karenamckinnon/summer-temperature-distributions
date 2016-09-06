% Main wrapper code for McKinnon et al (2016)
% The changing shape of Northern Hemisphere summer temperature distributions
% JGR Atmospheres
% doi: 10.1002/2016JD025292
  
figure('visible','off')

% user needs to fill these in
ghcnd_matpath = ''
dataDir = ''
cacheDir = '';
figDir = '';  
stationList = ''; % list of ISD stations produced from getHourly.R

yearStart = 1980;
yearEnd = 2015; 
months = [7 8]; 
percentiles = [5:5:95]; % percentiles for QR calculation
varname = 'TMIN' % 'TMAX' or 'TMIN'
% if runNumber = 0, components of calculations have already been done, and will be combined
% if runNumber > 0, then code wants to be parallelized (lots of computation) 
runNumber = 0; 
njitter = 100; % number of times to jitter data within precision bounds
nboot = 1000; % number of times to bootstrap to assess significance
nIntervals = 128; % how many processors do you want to use? 
latRange = [30 80]; % domain for analysis
lonRange = [-180 360];
sigCutoff = 0.95; % significance level for plotting (against null distribution assuming stationarity)

%% Fig. 1 (maps of all stations that pass QC globally) 
% show skewness and excess kurtosis
mapStations('TMAX', ghcnd_matpath, yearStart, yearEnd, months, [-90 90], ...
	[-180 360], cacheDir, figDir)

mapStations('TMIN', ghcnd_matpath, yearStart, yearEnd, months, [-90 90], ...
	[-180 360], cacheDir, figDir)
	

%% Get GHCND data, as well as QR results
[GHCND QRghcnd] = mainGHCND(yearStart, yearEnd, months, percentiles, varname, ...
	ghcnd_matpath, latRange, lonRange, ...
	runNumber, njitter, nboot, nIntervals);
	 
% Fig. 2 (Plot climatological seasonal cycle across stations)
plotSeasonalCycle(GHCND, figDir, varname)
	
%% Get null distribution based on stationarity assumption
QRnull = QRwrapperNull(GHCND, varname, nIntervals, months, runNumber, percentiles, ...
	cacheDir, nboot);
	
% Use the 'legendre' basis for results presented in paper 
basisName = 'legendre';
makePlot = 0; % plot bases? 
P = getBasis(percentiles, [], basisName, makePlot, figDir, []);

% if using PCA to get bases, need to input relevant QR and varname
% P = getBasis(percentiles, QR, basisName, makePlot, figDir, varname)
	
%% Project null distribution onto basis functions
[PCsNullRange, PCsNull] = bootstrapSig(GHCND, P, QRnull, varname, months, percentiles, ...
	cacheDir, nboot, sigCutoff);

%% Figs. 4-7 (trends in the 5, 50, and 95 percentiles)
datasource = 'GHCND';
% Make plots of changes in different percentiles
latRanges = [25 60; 25 75];% 25 72; 30 80];
lonRanges = [-130 -60; -15 180];% 80 180; -180 180];
cRanges = [20 50;-1.5 1.5];
regions = char('US','Eurasia'); 
percentilesToPlot = [5 50 95];

type = 'trend';
for ii = 1:2
	for jj = 1:length(percentilesToPlot)
		nameTag = ['GHCND.' varname '.month.' num2str(months(1)) '.' strtrim(regions(ii,:)) '.' num2str(percentilesToPlot(jj))];
		plotPercentiles(percentilesToPlot(jj), QRghcnd, GHCND, sigCutoff, varname, months, ...
			percentiles, figDir, nameTag, latRanges(ii, :), lonRanges(ii, :), cRanges, type)
	end  
end

%% Figs. 8, 10-13 (projections onto each basis)
type = 'slope';
for ii = 1:2
	nameTag = ['GHCND.' strtrim(regions(ii,:))]
	legendreBasis(QRghcnd, GHCND, P, percentiles, PCsNullRange, figDir,...
		datasource, varname, months, sigCutoff, percentilesToPlot,...
		latRanges(ii, :), lonRanges(ii, :), nameTag)
end

	
%%%%% HOURLY ANALYSIS %%%%%%%%%
	
%% Get hourly data from ISD for comparison
% need to change a few options
timesToSample = [7 17];
timeForQR = 7; % only perform QR on one time of day
nIntervalsISD = 64;
latRanges = [25 60; 25 72; 25 72];
lonRanges = [-130 -60; -10 80; 80 180];

[ISD QRisd] = mainHourly(yearStart, yearEnd, months, percentiles, varname, stationList, ...
	timesToSample, timeForQR, runNumber, nIntervalsISD, njitter, nboot);


%% Fig. S2 and S3 (plots comparing the results from GHCND to those of ISD)
compareIsdGhcndProj(QRghcnd, QRisd, GHCND, ISD, P, latRanges(1:3, :), lonRanges(1:3, :), varname, yearStart, yearEnd, figDir)

%% Fig. 14 (shifted distributions)
if strcmp(varname, 'TMAX')
	lats = [35 50; 45 60];
	lons = [-130 -115; 30 45];
	names = {'W.US','E.Europe'};
else
	lats = [38 48];
	lons = [-90 -70];
	names = {'NE.US'};
end
exampleDistributions(GHCND, lats, lons, names, varname, figDir, yearStart, yearEnd)





