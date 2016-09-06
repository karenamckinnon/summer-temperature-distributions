function QR = doQRnull(GHCND, months, percentiles, varname, nboot, timeSubset, Tsubset, cacheDir)
% perform QR after bootstrapping to calculate results given a null hypothesis of a stationary distribution

disp(['Calculating null QR for ' varname ''])

nyrs = length(unique(year(timeSubset)));
ndays = length(timeSubset)/nyrs;
nstations = size(Tsubset, 2);

years = year(timeSubset);
yearsCentered = (years - mean(years))';

QR.percentiles = percentiles;
OPTS.varname = varname;
OPTS.months = months;
QR.opts = OPTS;

betaNull = NaN(nstations, length(percentiles), nboot);

for ct = 1:nstations

	T0 = Tsubset(:, ct);
	% check to see if this station has already been calculated
	savename = getHash([cacheDir '/QRstation'], varname,...
		nboot, percentiles, months, T0);
	if exist(savename, 'file')
		disp(['Loading station ' num2str(ct) ' of ' num2str(nstations)])
		load(savename)
	else

		tic
		disp(['Calculating station ' num2str(ct) ' of ' num2str(nstations)])

		predictors0 = [ones(size(T0)), yearsCentered]; % num obs x 2

		% calculate quantization level
		if isfield(GHCND, 'id')
			T0 = 10*T0;
			[qlevs] = estimate_quantization_v2(timeSubset, T0); % runs on units of 1/10 C
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
		elseif isfield(GHCND, 'USAF') % ISD, assume at 0.1C precision based on data
			ci = [-0.05*(ones(length(T0), 1)) 0.05*(ones(length(T0), 1))];
		else % GHCND, all at 1F
			ci = 5/9*[-0.5*(ones(length(T0), 1)) 0.5*(ones(length(T0), 1))];
		end

		dQ0 = ci(:, 2) - ci(:, 1);
		mQ0 = mean(ci, 2);
 
		Tjitter = T0 + dQ0.*rand(size(dQ0)) - dQ0/2 + mQ0;

		% randomize data, keeping years together
		betaNullTMP = NaN(length(percentiles), nboot);
		for kk = 1:nboot

			index = randi(nyrs,1,nyrs);
			for pct = 1:numel(percentiles)

				% put back into matrix form for block boostrapping
				predmat = reshape(yearsCentered,[ndays nyrs]);
				valmat = reshape(Tjitter,[ndays nyrs]);

				% randomize by years
				V = valmat(:, index);

				ind = ~isnan(V);
				if length(find(ind)) > 1 % if not, just continue through bootstrap loop

					V = V(ind);
					P = [ones(size(V)), predmat(ind)];

					dum = rq(P, V, percentiles(pct)/100);
					betaNullTMP(pct, kk) = dum(2);
				end

			end


		end

		timestamp = datestr(now);
		save(savename, 'betaNullTMP','timestamp')
		toc

	end
	betaNull(ct, :, :) = betaNullTMP;

end

QR.betaNull = betaNull;







