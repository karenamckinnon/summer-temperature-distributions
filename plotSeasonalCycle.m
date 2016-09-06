function plotSeasonalCycle(GHCND, figDir, varname)

% calculate climatology 
doy = GHCND.time - datenum(year(GHCND.time), 1, 1) + 1;
 
% find leap days, and remove
notLeap = doy ~= 366;
T = GHCND.(varname)(notLeap, :);
timeNoLeap = GHCND.time(notLeap);
doyNoLeap = doy(notLeap);
 
% pull out complete years
idxStart = find(doyNoLeap == 1, 1, 'first');
idxEnd = find(doyNoLeap == 365, 1, 'last');
T = T(idxStart:idxEnd,:);
nyrs = size(T, 1) / 365;
Tmat = reshape(T, [365 nyrs size(T, 2)]);
% remove the mean from each year
TmatDev = bsxfun(@minus, Tmat, nanmean(Tmat, 1));
% and average 
TmatAvg = squeeze(nanmean(TmatDev, 2));
nstations = size(TmatAvg, 2);
 
tBasis = ([1:365]' - 0.5)/365;
basis1 = exp(2*pi*i*tBasis);
basis2 = exp(4*pi*i*tBasis);
basis3 = exp(6*pi*i*tBasis);
 
% project onto bases
eof1 = 2/length(tBasis)*sum(bsxfun(@times, TmatAvg, basis1), 1);
eof2 = 2/length(tBasis)*sum(bsxfun(@times, TmatAvg, basis2), 1);
eof3 = 2/length(tBasis)*sum(bsxfun(@times, TmatAvg, basis3), 1);

rec = real(bsxfun(@times, conj(eof1), basis1)) + ...
	real(bsxfun(@times, conj(eof2), basis2)) + ...
	real(bsxfun(@times, conj(eof3), basis3));
	
anomFromClim = TmatDev - reshape(repmat(rec,[nyrs 1]), [365 nyrs nstations]);
interannualSpread = prctile(anomFromClim, [2.5 97.5], 2);
	
% normalize everything by its own standard deviation
recNorm = bsxfun(@times, rec, 1./std(rec, [], 1));
spreadNorm1 = bsxfun(@times, squeeze(interannualSpread(:, 1, :)), 1./std(rec, [], 1));
spreadNorm2 = bsxfun(@times, squeeze(interannualSpread(:, 2, :)), 1./std(rec, [], 1));

recRange = prctile(recNorm, [0.5 99.5], 2)';
x = 1:365;
clf
% h = patch([x fliplr(x)], [recRange(1, :) fliplr(recRange(2, :))], 0.9*[1 1 1]);
h = patch([x fliplr(x)], [(mean(recNorm, 2) + mean(spreadNorm1, 2))' fliplr((mean(recNorm, 2) + mean(spreadNorm2, 2))')], 0.9*[1 1 1]);
set(h,'edgecolor','none')
% alpha(0.5)
hold on
plot(x, mean(recNorm, 2), 'k' ,'linewidth',2)
june1 = 152;
july1 = 182;
august31 = 243;
plot([june1 june1],[-4 3],'-.k')
plot([july1 july1],[-4 3],'--k', 'linewidth', 2)
plot([august31 august31],[-4 3],'--k', 'linewidth', 2)
xlim([1 365])
ylim([-3.5 2.5])
set(gca,'box','on','linewidth',1)
set(gca,'fontsize',10)
set(gca,'layer','top')

set(gcf, 'Color', 'w');
figName = [figDir '/averageSeasonalCycle.Spread.' varname '.pdf'];
export_fig(figName, '-m5','-a1')



