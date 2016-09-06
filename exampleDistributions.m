function exampleDistributions(GHCND, lats, lons, names, varname, figDir, yearStart, yearEnd)
% for a small subset of regions, plot the before and after distributions
   
for ct = 1:size(lats, 1)
 
	figName = ['shiftedDistr.' varname '.' names{ct} '.png'];

	locIdx = GHCND.loc(:,1)>lons(ct,1) & GHCND.loc(:,1)<lons(ct,2) & ...
		GHCND.loc(:,2)>lats(ct,1) & GHCND.loc(:,2)<lats(ct,2);
		
	disp([names{ct}])
	disp(['Number of stations going into plot: ' num2str(sum(locIdx)) ''])

	%{
	timeIdx = ismember(month(GHCND.time),7:8) & ismember(year(GHCND.time),yearStart:yearEnd);

	vals = GHCND.(varname)(timeIdx, locIdx);

	timeBlockSize = 62*10;

	vals1 = vals(1:timeBlockSize, :);
	vals2 = vals((end-timeBlockSize+1):end,:);

	% get mean of first period and remove
	meanVals1 = nanmean(vals, 1);
	anom1 = bsxfun(@minus, vals1, meanVals1);
	anom2 = bsxfun(@minus, vals2, meanVals1);

	minVal = min(min(vals1(:)), min(vals2(:)));
	maxVal = max(max(vals1(:)), max(vals2(:)));

	xi = linspace(-15,15, 101);
	xi = [xi(1:find(xi == -5.1)) -5 xi(find(xi == -5.1) + 1:end)];
	xi = [xi(1:68) 5 xi(69:end)];
	% add in 5's

	[y1, ~, bw1] = ksdensity(anom1(:), xi, 'bandwidth', 0.5);
	[y2, ~, bw2] = ksdensity(anom2(:), xi, 'bandwidth', 0.5);

	clf;
	% add shading for extremes (cold / warm)
	% thresholds: 5, 7.5
	xCold = find(xi == -7.5);
	Cold = find(xi == -5);
	Hot = find(xi == 5);
	xHot = find(xi == 7.5);
	colorxCold = rgb('MediumBlue');
	colorCold = rgb('CornflowerBlue');
	colorHot = rgb('LightCoral');
	colorxHot = rgb('FireBrick');
  
	h = patch([xi(1:Cold) fliplr(xi(1:Cold))], ...
		[zeros(size(xi(1:Cold))) fliplr(y2(1:Cold))], ...
		colorCold);
	set(h,'edgecolor','none');
	hold on
	h = patch([xi(1:xCold) fliplr(xi(1:xCold))], ...
		[zeros(size(xi(1:xCold))) fliplr(y2(1:xCold))], ...
		colorxCold);
	set(h,'edgecolor','none');
	h = patch([xi(Hot:end) fliplr(xi(Hot:end))], ...
		[zeros(size(xi(Hot:end))) fliplr(y2(Hot:end))], ...
		colorHot);
	set(h,'edgecolor','none');
	h = patch([xi(xHot:end) fliplr(xi(xHot:end))], ...
		[zeros(size(xi(xHot:end))) fliplr(y2(xHot:end))], ...
		colorxHot);
	set(h,'edgecolor','none');

	% get maxima
	[~, ind1] = max(y1);
	[~, ind2] = max(y2);

	plot([xi(ind1) xi(ind1)],[0 0.15],'--','color',0.7*[1 1 1], 'linewidth', 2);
	plot([xi(ind2) xi(ind2)],[0 0.15],'--k', 'linewidth', 2);

	plot(xi,y1,'color',0.7*[1 1 1],'linewidth',2)
	hold on
	plot(xi,y2,'k','linewidth',2)

	set(gca,'fontsize',12)

	orient landscape
	set(gca,'box','on')
	set(gca,'layer','top')
	grid on
	set(gcf,'color','w')
	ax = gca;
	ax.Position(4) = 0.4;
	ax.Position(3) = 0.775;
	ylim([0 0.12])
	xlim([-15 15])
	export_fig([figDir '/' figName], '-m3','-a1')
	%}

end