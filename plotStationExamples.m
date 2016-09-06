% Makes Fig. 3 -- example QR trends

sig = QRghcnd.fracSig > 0.95;
betaSig = QRghcnd.beta(:, :, 2);
betaSig(~sig) = NaN;
  
% calculate average summer skewness
timeInd = ismember(month(GHCND.time), months);
yrs = yearStart:yearEnd;
nyrs = length(yrs);
ndays = sum(timeInd)/nyrs;
nstations = size(GHCND.loc, 1);
% reshape into days x yrs x stations
Tmat = reshape(GHCND.(varname)(timeInd, :), [ndays nyrs nstations]);
skewEst = squeeze(nanmean(skewness(Tmat, [], 1), 2));
timeSubset = GHCND.time(timeInd);

stations1 = betaSig(:, percentiles == 5) > 0 & ...
	betaSig(:, percentiles == 50) > 0 & ...
	betaSig(:, percentiles == 95) > 0 & ...
	betaSig(:, percentiles == 5) > betaSig(:, percentiles == 50) & ...
	betaSig(:, percentiles == 95) > betaSig(:, percentiles == 5) & ...
	skewEst < 0;
	
% use a station that is in North America (idx 51)
stations1 = find(stations1);
stations1 = stations1(51);
	

% (2)
stations2 = betaSig(:, percentiles == 5) > 0 & ...
	betaSig(:, percentiles == 50) < 0 & ...
	betaSig(:, percentiles == 95) < 0 & ...
	skewEst > 0;
	
stations2 = find(stations2);
	
	
% plot
Xcentered = [yrs - mean(yrs)]';
percentilesToPlot = [5 50 95];
for ct = 1:2
	if ct == 1, idx = stations1; else, idx = stations2; end
	TmatUse = Tmat(:, :, idx);
	
	% jitter once for display
	T0 = 10*TmatUse(:); 
	[qlevs] = estimate_quantization_v2(timeSubset, T0);
	if size(qlevs, 1) == 1, qlevs = qlevs'; end
	ci = [-qlevs/2 qlevs/2];
	loc = find(qlevs==50/9);
	[a b] = dblround_ci( T0(loc) );
	nanind = isnan(a+b); % check to make sure no NaNs
	loc(nanind) = [];a(nanind) = [];b(nanind) = [];

	ci(loc, 1) = a-T0(loc);
	ci(loc, 2) = b-T0(loc);

	ci = ci/10; % switch back to C
	T0 = T0/10;
	
	dQ = ci(:, 2) - ci(:, 1);
	mQ = mean(ci, 2);
	T = T0 + dQ.*rand(size(dQ)) - dQ/2 + mQ;

	T = reshape(T,size(TmatUse));
	
	yfit = squeeze(QRghcnd.beta(idx, :, :))*[ones(size(Xcentered)) Xcentered]';
	clf
	plot(yrs, TmatUse, 'ok')
	hold on
	plot(yrs, T, '.k')
	plot(yrs, yfit(ismember(percentiles, 50), :),'color','k','linewidth', 2);
	%plot(yrs, yfit(ismember(percentiles, 50), :),'color','w','linewidth', 1);
	
	plot(yrs, yfit(ismember(percentiles, [5 95]), :),'-.','color','k','linewidth', 2);
	%plot(yrs, yfit(ismember(percentiles, [5 95]), :),'--','color','w','linewidth', 1);
	
	%plot(yrs, yfit(ismember(percentiles, 50), :),'color',0.5*[1 1 1],'linewidth', 2);
	%plot(yrs, yfit(ismember(percentiles, [5 95]), :),'--','color',0.5*[1 1 1],'linewidth', 2);
 
	xlim([1980 2015])
	
	set(gca,'fontsize',10)
	orient landscape
	
	set(gcf, 'Color', 'w');

	figName = [figDir '/exampleTS.' GHCND.id(idx,:) '.pdf'];
	export_fig(figName, '-m5','-a1')
	
	% make a map of the locatioin
	clf
	m_proj('miller', 'lat', [25 60], 'lon', [-130 -60]);
	m_coast('patch', 0.75*[1 1 1],'edgecolor','none');
	m_grid('box','on','linewidth', 2, 'linestyle','none','fontsize',10,'tickdir','in','xticklabels',[],'yticklabels',[]);
	hold on
	[sx sy] = m_ll2xy(GHCND.loc(idx,1), GHCND.loc(idx,2));
	scatter(sx,sy,150,'k','p','filled')
	set(gca,'fontsize', 10)

	orient landscape
	set(gcf, 'Color', 'w');

	figName = [figDir '/exampleTSmap.' GHCND.id(idx,:) '.pdf'];
	export_fig(figName, '-m5','-a1')

end

