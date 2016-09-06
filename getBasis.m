function P = getBasis(percentiles, QR, basisName, makePlot, figDir, varname)
  
% first, calculate basis for specified shifts in the moments
% this gives the normalization for all other bases

x = percentiles'/1e2;
k = 1e6;

%% Bases from specified shifts in a distribution using the Pearson system
P0 = 1*ones(size(x)); 
P1 = norminv(x, 0, 1.5) - norminv(x, 0, 1);  % chance in variance of 0.5

% Get difference in empirical inverse cdfs
D21 = pearsrnd(0, 1, -0.5, 3, k, 1); % change in skewness of 1
D22 = pearsrnd(0, 1, 0.5, 3, k, 1);
P2 = prctile(D22, percentiles') - prctile(D21, percentiles');
D31 = pearsrnd(0, 1, 0, 2, k, 1); % change in kurtosis of 1
D32 = pearsrnd(0, 1, 0, 4, k, 1);
P3 = prctile(D32, percentiles') - prctile(D31, percentiles');

% variance and kurtosis functions are colinear, whereas all others are ~orthogonal
P = [P0 P1 P2 P3];
stdVals = std(P, [], 1);

% overwrite these bases if using another type of basis

if strcmp(basisName, 'legendre') % general basis: Legendre polynomials

	legX = 2*(percentiles'/1e2) - 1;
	n = numel(legX);
	P0 = ones(n, 1);
	P1 = legX;
	P2 = 1/2*(3*legX.^2 - 1);
	P3 = 1/2*(5*legX.^3 - 3*legX);

	P = [P0 P1 P2 P3];

	% normalize so that the PCs have the same variance as those calculated using 
	% specified shifts in the distributions above
	
	P(:, 2:4) = bsxfun(@times, P(:, 2:4), stdVals(2:end)./std(P(:, 2:4), [], 1));

elseif strcmp(basisName, 'PCA') % specific basis: the PCAs for the dataset

	M = 10*QR.beta(:, :, 2)'; % change per decade
	n = numel(x);
	P = ones(n,1);

	nstations = size(M, 2);

	b = NaN(nstations,1);
	for ct = 1:nstations
		b(ct) = regress(M(:, ct),P);
	end

	% get variance explained
	Mrec = P*b';
	varexp0 = xcPH(M(:),Mrec(:));

	% remove this signal from M
	MRes = M - Mrec;

	% Perform PCA on the rest
	[U, S, V] = svd(MRes);
	SCut = S(1:n, 1:n);
	varexp = diag(SCut.^2)/trace(SCut.^2);

	NpcUse = 3;

	P = [ones(size(U(:,1))) U(:, 1:NpcUse)];
	if P(percentiles == 5, 2) > P(percentiles == 95, 2), P(:, 2) = -P(:, 2); end
	if P(percentiles == 50, 3) > P(percentiles == 95, 3), P(:, 3) = -P(:, 3); end
	if P(percentiles == 70, 4) > P(percentiles == 95, 4), P(:, 4) = -P(:, 4); end
 
	% normalize so that the PCs have the same variance as those calculated using 
	% specified shifts in the distributions above
	
	P(:, 2:4) = bsxfun(@times, P(:, 2:4), stdVals(2:end)./std(P(:, 2:4), [], 1));
end

if makePlot
	clf
	hold on
	or = get(gca,'colororder');
	or(3,:) = [];
	set(gca,'colororder',or)
	for ct = 1:size(P,2)
		plot(percentiles, P(:,ct), 'linewidth', 2);
	end
	ylim([-1.25 1.25])
	set(gca,'fontsize',12)

	if strcmp(basisName, 'legendre')
		hleg = legend('P0','P1','P2','P3','location','southeast');
	elseif strcmp(basisName, 'PCA')
		hleg = legend('B0','B1','B2','B3','location','southeast');
	elseif strcmp(basisName, 'moments')
		hleg = legend('mean','var','skew','kurt','location','southeast');
	end
	set(hleg,'fontsize',12)

	orient landscape
	ax = gca;
	ax.Position(4) = 0.6;
	ax.Position(3) = 0.775;
	grid on
	orient landscape
	set(gcf,'color','w')
	set(gca,'box','on')
	set(gca,'layer','top')
	if strcmp(basisName, 'PCA')
		export_fig([figDir '/basis.' basisName '.' varname '.png'],'-m3','-a1')
	else
		export_fig([figDir '/basis.' basisName '.png'],'-m3','-a1')
	end


end