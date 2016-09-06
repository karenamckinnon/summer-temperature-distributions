function compareIsdGhcndProj(QRghcnd, QRisd, GHCND, ISD, P, latRanges, lonRanges, varname, yearStart, yearEnd, figDir)
% Compare the projections of the QR trends from ghcnd and isd to each other
  
% find distance between points
deg2rad = pi/180; 
RR=6378.137; % radius of the earth in km
dlon = deg2rad*bsxfun(@minus, ISD.loc(:, 1), GHCND.loc(:, 1)');
dlat = deg2rad*bsxfun(@minus, ISD.loc(:, 2), GHCND.loc(:, 2)');
a = (sin(dlat/2)).^2 + bsxfun(@times, cos(deg2rad*GHCND.loc(:, 2))', cos(deg2rad*ISD.loc(:, 2))) .* (sin(dlon/2)).^2;
c = 2 * atan2( sqrt(a), sqrt(1-a) );
d = RR * c;

% for each ISD station, find min distance to a GHCND station
[minDist, distInd] = min(d, [], 2);

% require that stations are within 100km of each other and that elevations are similar
cutoff = 100;
idx = find(minDist < cutoff);
distInd = distInd(idx);
minDist = minDist(idx);

% get paired locations
ghcndLocs = GHCND.loc(distInd, :);
isdLocs = ISD.loc(idx, :);
  
% remove locations with dz > 100 m
dz = ghcndLocs(:, 3) - isdLocs(:, 3);
cutoffZ = 100;
idxZ = abs(dz) < cutoffZ;
ghcndLocs = ghcndLocs(idxZ, :);
isdLocs = isdLocs(idxZ, :);

% get original indices
isdIdx = idx(idxZ);
ghcndIdx = distInd(idxZ);

% compare QR values for all percentiles
betaISD = QRisd.beta(isdIdx, :, :);
betaGHCND = QRghcnd.beta(ghcndIdx, :, :);

% if stations have differing mean values, discard
deltaMean = mean(sqrt((betaISD(:, :, 1) - betaGHCND(:, :, 1)).^2), 2);

cutoffMean = 1;
idxMean = find(deltaMean < cutoffMean);
  
isdIdx = isdIdx(idxMean);
ghcndIdx = ghcndIdx(idxMean);

% compare QR values for all percentiles
betaISD = 10*QRisd.beta(isdIdx, :, 2)';
betaGHCND = 10*QRghcnd.beta(ghcndIdx, :, 2)';

nstations = size(betaISD, 2)


% project these onto Legendre basis functions
PCsISD = NaN(nstations, size(P, 2));
PCsGHCND = PCsISD;
for ct = 1:nstations
	PCsISD(ct,:) = regress(betaISD(:,ct),P); % nstation x PC
	PCsGHCND(ct, :) = regress(betaGHCND(:, ct), P);
end

% Make a map of where the stations are

clf
m_proj('miller','lon',[-180 180],'lat',[20 80]);
m_coast('patch', 0.75*[1 1 1],'edgecolor','none');
m_grid('box','on','linewidth', 0.5, 'linestyle','none','fontsize',8,'tickdir','in');
hold on
[sx1 sy1] = m_ll2xy(GHCND.loc(ghcndIdx, 1), GHCND.loc(ghcndIdx, 2));
[sx2 sy2] = m_ll2xy(ISD.loc(isdIdx, 1), ISD.loc(isdIdx, 2));
scatter(sx1, sy1, 3, 'k', 'o', 'filled')
scatter(sx2, sy2, 2, 0.9*[1 1 1], 'p','filled')
set(gca,'fontsize',10)
 
set(gcf, 'Color', 'w');
figName = [figDir '/compare-ghcnd-isd-' varname '-' num2str(yearStart) '-' num2str(yearEnd) '.png'];
export_fig(figName, '-m4','-a1')
 
 
%% Compare the QR values, and use different regions
nRegions = size(latRanges, 1);
nBasis = size(P,2);
clf 
for ii = 1:nRegions
	regionIdx = GHCND.loc(ghcndIdx, 1) > lonRanges(ii, 1) & ...
		GHCND.loc(ghcndIdx, 1) < lonRanges(ii, 2) & ...
		GHCND.loc(ghcndIdx, 2) > latRanges(ii, 1) & ...
		GHCND.loc(ghcndIdx, 2) < latRanges(ii, 2);
	for jj = 1:nBasis
		plotIdx = (jj - 1)*nRegions + ii;
		% subplot(nBasis, nRegions, plotIdx);
		clf
		
		plot([-1.5 1.5],[-1.5 1.5],'--','color',0.7*[1 1 1], 'linewidth', 2)
		hold on
		scatter(PCsGHCND(regionIdx, jj), ...
			PCsISD(regionIdx, jj), 25, 'k','filled','o')
		
		xlim([-1.5 1.5])
		ylim([-1.5 1.5])
		axis square
		
		% calculate bias (Celsius per decade units)
		b = mean(PCsGHCND(regionIdx, jj)) - mean(PCsISD(regionIdx, jj));
		text(-1.4, 1.35, [num2str(1/100*round(100*b))],'fontsize',20)
		set(gca,'linewidth',1)
		
		grid on
		
		set(gca,'xtick',-1:1:1)
		set(gca,'xticklabel',-1:1:1);
		set(gca,'ytick',-1:1:1)
		set(gca,'yticklabel',-1:1:1);
		set(gca,'fontsize',16)
		set(gcf, 'Color', 'w');
		figName = [figDir '/ghcnd.vs.isd.projections.region.' num2str(ii) '.basis.' num2str(jj) '.' varname '.' num2str(yearStart) '.' num2str(yearEnd) '.png'];
		export_fig(figName, '-m2','-a1')
	
	end
end
	




