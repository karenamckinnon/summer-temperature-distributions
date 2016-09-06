function [PCsNullRange, PCsNull] = bootstrapSig(GHCND, P, QRnull, varname, months, percentiles, ...
	cacheDir, nboot, sigCutoff)
% get quantile regression trends based on null hypothesis that fields are stationary
  

% get distribution of legendre basis functions
savename = getHash([cacheDir '/PCsnull'], varname, months, percentiles, ...
	P, nboot, GHCND.id);

if exist(savename, 'file')
	load(savename)
else
	nstations = size(QRnull.betaNull, 1);
	PCsNull = NaN(nstations, size(P, 2), nboot);
	for kk = 1:nboot

		M = 10*QRnull.betaNull(:, :, kk)'; % change per decade

		% get variance explained by Legendre polynomials
		for ct = 1:nstations
			PCsNull(ct,:,kk) = regress(M(:,ct),P);
		end

	end
	timestamp = datestr(now);
	save(savename,'timestamp','PCsNull')
end

% calculate range that is consistent with null hypothesis
% use a two-sided test
prcRanges = 100*[(1-sigCutoff)/2 (1+sigCutoff)/2];

PCsNullRange = mean(abs(prctile(PCsNull, prcRanges, 3)),3);






