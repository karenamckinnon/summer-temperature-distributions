function mapStations(varname, ghcnd_matpath, yearStart, yearEnd, months, latRange, ...
	lonRange, cacheDir, figDir)
% make a map of stations that pass QC
 
GHCNDdummy = getGHCND(varname, ghcnd_matpath, yearStart, yearEnd, months, latRange, lonRange, cacheDir);
GHCNDdummy = doQC(GHCNDdummy, months, varname);
 
% calculate average summer skewness
timeInd = ismember(month(GHCNDdummy.time), months);
nyrs = length(yearStart:yearEnd);
ndays = sum(timeInd)/nyrs;
nstations = size(GHCNDdummy.loc, 1);
% reshape into days x yrs x stations
Tmat = reshape(GHCNDdummy.(varname)(timeInd, :), [ndays nyrs nstations]);
skewEst = squeeze(nanmean(skewness(Tmat, [], 1), 2));
kurtEst = squeeze(nanmean(kurtosis(Tmat, [], 1), 2));

% switch to excess kurtosis
kurtEst = kurtEst - 3;

% test for normality via bootstrapping what range of skewness could be expected with a normal distribution

nboot = 1e5;
sample = randn(ndays, nyrs, nboot);
skewEstNormal = squeeze(mean(skewness(sample,[],1),2));
kurtEstNormal = squeeze(mean(kurtosis(sample,[],1),2)) - 3;

%  95% range of possible skewnesses
skewCIs = prctile(skewEstNormal, [2.5 97.5])
kurtCIs = prctile(kurtEstNormal, [2.5 97.5])


cmap = flipud(lbmap(20, 'BrownBlue'));
cmap(10:11, :) = 0.9*[1 1 1; 1 1 1];

% plot stations
clf
 
m_proj('miller','lon',[-180 180],'lat',[-60 80]);
m_coast('patch', 0.7*[1 1 1],'edgecolor','none');
m_grid('box','on','linewidth', 2, 'linestyle','none','fontsize',10,'tickdir','in');
hold on
[sx sy] = m_ll2xy(GHCNDdummy.loc(:, 1), GHCNDdummy.loc(:, 2));
scatter(sx, sy, 2, skewEst, 'o', 'filled')
% plot lines for 30 to 72

hold on
m_plot([-180 180], [30 30], 'k','linewidth',1)
m_plot([-180 180], [72 72], 'k','linewidth',1)
h = colorbar;
ylabel(h, 'Skewness','fontsize',10)
colormap(cmap)
caxis([-1 1])
set(gca,'fontsize',10)
set(h,'fontsize',10)

set(gcf, 'Color', 'w');
figName = [figDir '/globalMap' varname '-' num2str(yearStart) '-' num2str(yearEnd) '-skewness.pdf'];
export_fig(figName, '-m5','-a1')

% span -0.5 to 1.5 in steps of 0.1
cmap = flipud(lbmap(20, 'BrownBlue'));
cmap(8:11, :) = repmat(0.9*[1 1 1],[4 1]);

% plot stations
clf

m_proj('miller','lon',[-180 180],'lat',[-60 80]);
m_coast('patch', 0.7*[1 1 1],'edgecolor','none');
m_grid('box','on','linewidth', 2, 'linestyle','none','fontsize',10,'tickdir','in');
hold on
[sx sy] = m_ll2xy(GHCNDdummy.loc(:, 1), GHCNDdummy.loc(:, 2));
scatter(sx, sy, 2, kurtEst, 'o', 'filled')
h = colorbar;
ylabel(h, 'Excess kurtosis','fontsize',10)
colormap(cmap)
caxis([-1 1])
set(gca,'fontsize',10)
set(h,'fontsize',10)

set(gcf, 'Color', 'w');
figName = [figDir '/globalMap' varname '-' num2str(yearStart) '-' num2str(yearEnd) '-excess-kurtosis.pdf'];
export_fig(figName, '-m5','-a1')