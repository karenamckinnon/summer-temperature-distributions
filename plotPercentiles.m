function plotPercentiles(percentileToPlot, QR, GHCND, sigCutoff, varname, months, ...
	percentiles, figDir, nameTag, latRange, lonRange, cRanges, type)
% Plot trends in specified percentiles
 
clf  
nclmaps
cmap = thermal_soft;
 
if strcmp(type,'delta')
	beta = QR.dbeta;
	fracSig = QR.dfracSig;
else 
	beta = QR.beta(:, :, 2);
	fracSig = QR.fracSig;
end

if any(~isnan(fracSig(:, percentiles == percentileToPlot)))
	sigStations = fracSig(:, percentiles == percentileToPlot) > sigCutoff;
else % if we haven't calculated significance, just plot everything as significant
	sigStations = isnan(fracSig(:, percentiles == percentileToPlot));
end

 
plotIdx = GHCND.loc(:,1) > min(lonRange) & GHCND.loc(:,1) < max(lonRange) & ...
	GHCND.loc(:,2) > min(latRange) & GHCND.loc(:,2) < max(latRange);
	
% how many stations, and how many sig?
disp(nameTag)
disp(['Total stations: ' num2str(sum(plotIdx)) ''])
disp(['Sig stations: ' num2str(sum(sigStations & plotIdx)) ''])


clf
m_proj('miller', 'lat', latRange, 'lon', lonRange);
m_coast('patch', 0.75*[1 1 1],'edgecolor','none');
m_grid('box','on','linewidth', 2, 'linestyle','none','fontsize',10,'tickdir','in');
 
hold on
[sx sy] = m_ll2xy(GHCND.loc(:,1), GHCND.loc(:,2)); 
scatter(sx(~sigStations & plotIdx), sy(~sigStations & plotIdx), 10, 10*beta(~sigStations & plotIdx, percentiles == percentileToPlot), 'filled', 'o', 'linewidth', 2);
scatter(sx(sigStations & plotIdx), sy(sigStations & plotIdx), 20, 10*beta(sigStations & plotIdx, percentiles == percentileToPlot), '+', 'linewidth', 2);

h = colorbar;
ylabel(h,'\circC/decade')
colormap(cmap)
set(gca,'fontsize', 10)

caxis(cRanges(2, :))
set(gca,'fontsize',10)
set(h,'fontsize',10)
orient landscape
set(gcf, 'Color', 'w');

figName = [figDir '/' nameTag '.pdf'];
export_fig(figName, '-m4','-a1')

