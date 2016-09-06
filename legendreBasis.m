function legendreBasis(QR, GHCND, P, percentiles, PCsNullRange, figDir,...
	datasource, varname, months, sigCutoff, percentilesToPlot,...
	latRange, lonRange, nameTag)
% Plot projection coefficients onto Legendre basis functions
  
M = 10*QR.beta(:, :, 2)'; % change per decade
nstations = size(M, 2);

nclmaps
cmap = thermal_soft;

% find locations that are in specified region
locIdx = GHCND.loc(:,1) > min(lonRange) & GHCND.loc(:,1) < max(lonRange) & ...
		GHCND.loc(:,2) > min(latRange) & GHCND.loc(:,2) < max(latRange);

% get variance explained by Legendre polynomials
PCs = NaN(nstations, size(P, 2));
for ct = 1:nstations
	PCs(ct,:) = regress(M(:,ct),P); % nstation x PC
end

sigVals = abs(PCs) > PCsNullRange;


% Calculate the fraction of stations in each region with pos/neg trends
% as well as sig/not sig
% for Asia, separate Japan
if ~isempty(strfind(nameTag, 'Asia'))

	japanIdx = ismember(GHCND.id(:, 1:2), 'JA', 'rows');

	for kk = 1:2

		if kk == 1
			locIdxUse = locIdx & japanIdx;
		else
			locIdxUse = locIdx & ~japanIdx;
		end

		sigVals = abs(PCs) > PCsNullRange;

		% calculate the fraction of stations that are pos/neg, as well as sig pos/neg
		% for each region

		fracSig = sum(sigVals(locIdxUse, :), 1)/sum(locIdxUse);
		fracPos = sum(PCs(locIdxUse, :) > 0)/sum(locIdxUse);
		fracNeg = sum(PCs(locIdxUse, :) < 0)/sum(locIdxUse);
		sigValsPos = sum(sigVals(locIdxUse, :) & PCs(locIdxUse, :) > 0)/sum(locIdxUse);
		sigValsNeg = sum(sigVals(locIdxUse, :) & PCs(locIdxUse, :) < 0)/sum(locIdxUse);

		if kk == 1
			disp(['Japan'])
		else
			disp([nameTag])
		end

		disp(['Total number of stations: ' num2str(sum(locIdxUse)) ''])
		disp(['Frac pos: ' num2str(fracPos) ''])
		disp(['Frag neg: ' num2str(fracNeg) ''])
		disp(['Frac sig pos: ' num2str(sigValsPos) ''])
		disp(['Frac sig neg: ' num2str(sigValsNeg) ''])

	end

else

	sigVals = abs(PCs) > PCsNullRange;

	% calculate the fraction of stations that are pos/neg, as well as sig pos/neg
	% for each region

	fracSig = sum(sigVals(locIdx, :), 1)/sum(locIdx);
	fracPos = sum(PCs(locIdx, :) > 0)/sum(locIdx);
	fracNeg = sum(PCs(locIdx, :) < 0)/sum(locIdx);
	sigValsPos = sum(sigVals(locIdx, :) & PCs(locIdx, :) > 0)/sum(locIdx);
	sigValsNeg = sum(sigVals(locIdx, :) & PCs(locIdx, :) < 0)/sum(locIdx);

	disp([nameTag])
	disp(['Total number of stations: ' num2str(sum(locIdx)) ''])
	disp(['Frac pos: ' num2str(fracPos) ''])
	disp(['Frag neg: ' num2str(fracNeg) ''])
	disp(['Frac sig pos: ' num2str(sigValsPos) ''])
	disp(['Frac sig neg: ' num2str(sigValsNeg) ''])

end


% get variance explained overall
Mrec = P*PCs';
varexp = xcPH(M(:),Mrec(:));
disp(['Variance explained for ' varname ': ' num2str(varexp) ''])

% get fraction of change attributed to each basis function
for ct = 1:size(P,2);
	val = P(:,ct)*PCs(:,ct)';
	MrecPC = val;
	% how much variance explained?
	varexpPC(ct) = xcPH(MrecPC(:, locIdx), M(:, locIdx));
end
disp(['Variance explained by bases: ' num2str(varexpPC) ''])
disp(['Fractional variance explained by last three bases: ' num2str(varexpPC(2:end)/(1-varexpPC(1))) ''])


% any by stations, per basis function
varexpStation = NaN(size(PCs));
for jj = 1:nstations

	for ii = 1:size(P, 2)
		PCtemp(ii) = regress(M(:, jj), P(:, ii));
		Mrec = P(:, ii)*PCtemp(ii);
		varexpStation(jj, ii) = xcPH(M(:, jj), Mrec);
	end

end

x = percentiles'/1e2;

% Make plots
% (1) Basis functions (Fig. 8)
%colorOrder = [rgb('SeaGreen'); rgb('CornflowerBlue');rgb('Crimson');rgb('Goldenrod')];
clf
%set(gca, 'ColorOrder', colorOrder, 'NextPlot', 'replacechildren');

plot(percentiles, P(:, 1), '-k', 'linewidth', 2)
hold on
plot(percentiles, P(:, 2), '--k', 'linewidth', 2)
plot(percentiles, P(:, 3), ':k', 'linewidth', 2)
plot(percentiles, P(:, 4), '-.k', 'linewidth', 2)
ylim([-1.25 1.25])
set(gca,'fontsize',16)

hleg = legend('P0','P1','P2','P3','location','southeast');
set(hleg,'fontsize',16)
  
orient landscape
set(gcf,  'PaperPositionMode', 'manual');
set(gcf,  'PaperUnits', 'inches');
set(gcf,  'PaperPosition', [0.25 1 10 5]);
saveas(gcf, [figDir '/legendreBasis.pdf']);

% (2) coefficients for each basis function
% plot values as a function of space

for ct = 1:size(P,2)

	clf
	m_proj('miller', 'lat', latRange, 'lon', lonRange);
	m_coast('patch', 0.75*[1 1 1],'edgecolor','none');
	m_grid('box','on','linewidth', 2, 'linestyle','none','fontsize',10,'tickdir','in');
	hold on
	[sx sy] = m_ll2xy(GHCND.loc(:,1), GHCND.loc(:,2));
	% original vals: 12, 8, 12, 12
	if ct == 1
		scatter(sx(locIdx & sigVals(:, ct)), sy(locIdx & sigVals(:, ct)), 20, PCs(locIdx & sigVals(:, ct), ct), '+', 'linewidth', 1);
		scatter(sx(locIdx & ~sigVals(:, ct)), sy(locIdx & ~sigVals(:, ct)), 16, PCs(locIdx & ~sigVals(:, ct), ct), 'filled', 'o', 'linewidth', 2);

	else

		scatter(sx(locIdx & sigVals(:, ct)), sy(locIdx & sigVals(:, ct)), 20.*varexpStation(locIdx & sigVals(:, ct), ct), PCs(locIdx & sigVals(:, ct), ct), '+', 'linewidth', 1);
		scatter(sx(locIdx & ~sigVals(:, ct)), sy(locIdx & ~sigVals(:, ct)), 16.*varexpStation(locIdx & ~sigVals(:, ct), ct), PCs(locIdx & ~sigVals(:, ct), ct), 'filled', 'o', 'linewidth', 2);
		
	end
	h = colorbar;
	colormap(cmap)
	caxis([-1 1])

	set(gca,'fontsize',10)
	set(h,'fontsize',10)
	orient landscape
	set(gcf, 'Color', 'w');

	figName = [figDir '/map.legendre.P' num2str(ct-1) '.' nameTag '.' varname '.month.' num2str(months(1)) '.png'];
	export_fig(figName, '-m4','-a1')
	

end











