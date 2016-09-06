function plotTrendCorrs
% Makes Fig. 9

figure('visible','off')

% user must fill these in
ghcnd_matpath = '';
dataDir = '';
cacheDir = '';
figDir = '';

yearStart = 1980;
yearEnd = 2015;
months = [7 8];
percentiles = [5:5:95];
runNumber = 0; % if runNumber > 0, code is going to analyze part of data
njitter = 100;
nboot = 1000;
nIntervals = 128;
latRange = [30 80]; 
lonRange = [-180 360];


% delta percentiles
dpercentiles = [95 50; 5 50; 95 5];

%% Get GHCND part
[~, QRTx] = mainGHCND(yearStart, yearEnd, months, percentiles, 'TMAX', ...
	ghcnd_matpath, latRange, lonRange, ...
	runNumber, njitter, nboot, nIntervals);


 %% Get GHCND part
[~, QRTn] = mainGHCND(yearStart, yearEnd, months, percentiles, 'TMIN', ...
	ghcnd_matpath, latRange, lonRange, ...
	runNumber, njitter, nboot, nIntervals);


Cx = corrcoef(QRTx.beta(:, :, 2));
Cn = corrcoef(QRTn.beta(:, :, 2));

Cx = tril(Cx);
Cn = triu(Cn);

Cx(Cx == 0) = Cn(Cx == 0);

clf
pcolorPH(percentiles, percentiles, Cx)
colorbar
caxis([0.5 1])
colormap(flipud(lbmap(10,'redblue')))
hold on
plot([0 100],[0 100],'-k','linewidth',2)
set(gca,'box','on')
set(gca,'layer','top')
set(gcf,'renderer','painters')
set(gca,'xtick',[10:10:90])

set(gca,'fontsize',10)
orient landscape
set(gcf, 'Color', 'w');

figName = [figDir '/corrTrends.pdf'];
export_fig(figName, '-m4','-a1')
