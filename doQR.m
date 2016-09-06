function QR = doQR(GHCND, months, percentiles, njitter, varname, stationCountsforQR, ...
    Tsubset, timeSubset, nboot, cacheDir)
% perform quantile regression on each station, for specified months
% View data from a given year (within the month spans) as interchangeable

disp(['Calculating QR for ' varname ''])

nyrs = length(unique(year(timeSubset)));
ndays = length(timeSubset)/nyrs;
nstations = size(Tsubset, 2);

years = year(timeSubset);
yearsCentered = (years - mean(years))';

% order for correlations if we are looking at 5, 50, 95 (or any 3)
if length(percentiles) == 3 % don't want to do cross correlation across all percentiles
    order = perms(1:length(percentiles));
    order = unique(order(:,1:2),'rows');
    order(order(:,1)>order(:,2),:) = [];
end

QR.percentiles = percentiles;
OPTS.varname = varname;
OPTS.months = months;
OPTS.njitter = njitter;
OPTS.stationCountsforQR = stationCountsforQR;
if length(percentiles) == 3, OPTS.order = order; end
QR.opts = OPTS;

beta = NaN(nstations, length(percentiles), 2);
betaSD = beta;
dbeta = NaN(nstations, 3);
bootstrapSD = NaN(nstations, length(percentiles));
bootstrapSD_OLS = NaN(nstations, 1);
dbootstrapSD = NaN(nstations, 3);
fracSig = NaN(nstations, length(percentiles));
dfracSig = NaN(nstations, 3);
 
if length(percentiles) == 3, rho = NaN(nstations, size(order,1)); end
for ct = 1:nstations
	T0 = Tsubset(:, ct);
	% check to see if this station has already been calculated
	savename = getHash([cacheDir '/QRstation'], varname,...
		njitter, nboot, percentiles, months, T0);
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
			% Quantization estimate for rounded data
			% See Rhines et al (2015) Decoding the precision of historical temperature observations, QJRMS for more details
			[qlevs] = estimate_quantization_v2(timeSubset, T0); % runs on units of 1/10 C
			if size(qlevs, 1) == 1, qlevs = qlevs'; end
			ci = [-qlevs/2 qlevs/2];
			loc = find(qlevs==50/9); % Fahrenheit measurements
			[a b] = dblround_ci( T0(loc) );
			nanind = isnan(a+b); % check to make sure no NaNs
			loc(nanind) = [];a(nanind) = [];b(nanind) = [];

			ci(loc, 1) = a-T0(loc);
			ci(loc, 2) = b-T0(loc);

			ci = ci/10; % switch back to C
			T0 = T0/10;
		elseif isfield(GHCND, 'USAF') % ISD, assume at 0.1C precision based on data
			ci = [-0.05*(ones(length(T0), 1)) 0.05*(ones(length(T0), 1))];
		else % USHCN, all at 1F
			ci = 5/9*[-0.5*(ones(length(T0), 1)) 0.5*(ones(length(T0), 1))];
		end

		dQ0 = ci(:, 2) - ci(:, 1);
		mQ0 = mean(ci, 2);

		% remove missing vals
		pl = ~isnan(T0);
		T = T0(pl);
		predictors = predictors0(pl, :);

		dQ = dQ0(pl, :);
		mQ = mQ0(pl);

		QRslopes = NaN(njitter, length(percentiles), 2);

		for ii = 1:njitter
			% Tjitter = T + (-jitter_width + (2*jitter_width)*rand(size(T)));
			Tjitter = T + dQ.*rand(size(dQ)) - dQ/2 + mQ; % variable jitter width
			for jj = 1:length(percentiles)
				QRslopes(ii, jj, :) = rq(predictors, Tjitter, percentiles(jj)/100);
			end
		end 

		dQRslopes(:, 1) = QRslopes(:, percentiles == 95, 2) - QRslopes(:, percentiles == 50, 2);
		dQRslopes(:, 2) = QRslopes(:, percentiles ==  5, 2) - QRslopes(:, percentiles == 50, 2);
		dQRslopes(:, 3) = QRslopes(:, percentiles == 95, 2) - QRslopes(:, percentiles ==  5, 2);

		% get best estimate as median across jittered samples
		% record standard deviation as well

		betaTMP = squeeze(median(QRslopes, 1));
		betaSDTMP = squeeze(std(QRslopes, [], 1));
		dbetaTMP = squeeze(median(dQRslopes, 1));

		clear QRslopes dQRslopes

		beta_boot = NaN(length(percentiles), nboot);
		beta_ols_boot = NaN(1, nboot);

		% also do bootstrapping
		for kk = 1:nboot
			% now jittering every time, because I'd like an estimate of complete error
			Tjitter = T0 + dQ0.*rand(size(dQ0)) - dQ0/2 + mQ0;
			index = randi(nyrs,1,nyrs);
			for pct = 1:numel(percentiles)

				% put back into matrix form for block boostrapping
				predmat = reshape(yearsCentered,[ndays nyrs]);
				valmat = reshape(Tjitter,[ndays nyrs]);
				% residuals from best estimate of slope

				bestFit = predictors0*betaTMP(pct,:)';
				resmat = reshape(Tjitter - bestFit,[ndays nyrs]);

				% jitter the residuals, and add original slope estimate back in
				V = resmat(:, index);
				V = V(:) + bestFit;

				ind = ~isnan(V);
				if length(find(ind)) > 1 % if not, just continue through bootstrap loop

					V = V(ind);
					P = [ones(size(V)), predmat(ind)];

					dum = rq(P, V, percentiles(pct)/100);
					beta_boot(pct,kk) = dum(2);
					b = regress(V,P);
					beta_ols_boot(kk) = b(2);
				end

			end


		end
		% also look at 95-50, 5-50, 95-5
		dbeta_boot(1, :) = beta_boot(percentiles == 95, :) - beta_boot(percentiles == 50, :);
		dbeta_boot(2, :) = beta_boot(percentiles ==  5, :) - beta_boot(percentiles == 50, :);
		dbeta_boot(3, :) = beta_boot(percentiles == 95, :) - beta_boot(percentiles ==  5, :);

		if length(percentiles) == 3 && nboot > 0
			for mm = 1:size(order,1)
				rhoTMP(mm) = xcPH(beta_boot(order(mm, 1), :),beta_boot(order(mm,2), :),1);
			end
		else
			rhoTMP = [];
		end

		if nboot > 0
			% calculate significance of trends
			% significance defined as fraction of booted trends that are greater than zero if trend positive, and the converse
			for pct = 1:numel(percentiles)
				if betaTMP(pct,2) >= 0
					fracSigTMP(pct) = length(find(squeeze(beta_boot(pct,:))>0))/length(find(~isnan(squeeze(beta_boot(pct,:)))));
				elseif betaTMP(pct,2) < 0
					fracSigTMP(pct) = length(find(squeeze(beta_boot(pct,:))<0))/length(find(~isnan(squeeze(beta_boot(pct,:)))));
				end
			end

			for pct = 1:3
				if dbetaTMP(pct) >= 0
					dfracSigTMP(pct) = length(find(squeeze(dbeta_boot(pct,:))>0))/length(find(~isnan(squeeze(dbeta_boot(pct,:)))));
	 			elseif dbetaTMP(pct) < 0
					dfracSigTMP(pct) = length(find(squeeze(dbeta_boot(pct,:))<0))/length(find(~isnan(squeeze(dbeta_boot(pct,:)))));
				end 
			end 

			% save standard deviation of bootstrap samples
			bootstrapSDTMP = nanstd(beta_boot, [], 2);
			bootstrapSD_OLSTMP = nanstd(beta_ols_boot);
			dbootstrapSDTMP = nanstd(dbeta_boot, [], 2);
		end
		timestamp = datestr(now);
		save(savename, 'betaTMP', 'betaSDTMP', 'rhoTMP', 'fracSigTMP', 'bootstrapSDTMP',...
			'bootstrapSD_OLSTMP', 'dbetaTMP', 'dfracSigTMP', 'dbootstrapSDTMP','timestamp')
		toc

	end
	beta(ct, :, :) = betaTMP;
	betaSD(ct, :, :) = betaSDTMP;
	dbeta(ct, :) = dbetaTMP;
	if length(percentiles) == 3 && nboot > 0
		rho(ct, :) = rhoTMP;
	else
		rho = [];
	end
	fracSig(ct,:) = fracSigTMP;
	dfracSig(ct, :) = dfracSigTMP;
	bootstrapSD(ct,:) = bootstrapSDTMP;
	bootstrapSD_OLS(ct) = bootstrapSD_OLSTMP;
	dbootstrapSD(ct,:) = dbootstrapSDTMP;
	clear betaTMP betaSDTMP rhoTMP fracSigTMP bootstrapSDTMP bootstrapSD_OLSTMP timestamp dbootstrapSDTMP dfracSigTMP dbetaTMP

end

QR.beta = beta;
QR.betaSD = betaSD;
QR.rho = rho;
QR.bootstrapSD = bootstrapSD;
QR.bootstrapSD_OLS = bootstrapSD_OLS;
QR.fracSig = fracSig;
QR.dfracSig = dfracSig;
QR.dbootstrapSD = dbootstrapSD;
QR.dbeta = dbeta;






