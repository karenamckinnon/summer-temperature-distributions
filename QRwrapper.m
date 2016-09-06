function QR = QRwrapper(GHCND, varname, nIntervals, months, runNumber, percentiles, ...
	njitter, cacheDir, nboot);

% Get station count
nstations = size(GHCND.(varname), 2);
% divide stations among number of machines that are available
intervalsQR = round(linspace(1, nstations, nIntervals + 1));
if intervalsQR(end) ~= nstations; intervalsQR(end) = nstations; end
timeInd = find(month(GHCND.time) >= min(months) & month(GHCND.time) <= max(months));
Tsubset0 = GHCND.(varname)(timeInd, :);

if runNumber > 0 % Calculate QR
  
	stationCountsforQR = [intervalsQR(runNumber): intervalsQR(runNumber + 1)];
	% do QR
	timeInd = find(month(GHCND.time) >= min(months) & month(GHCND.time) <= max(months));
	Tsubset = Tsubset0(:, ismember(1:nstations, stationCountsforQR));
	timeSubset = GHCND.time(timeInd);

	savename = getHash([cacheDir '/QRstationSet'], varname, months, percentiles, njitter, Tsubset, ...
		stationCountsforQR, nboot);

	if exist(savename, 'file')
		load(savename)
	else
		QRpartial = doQR(GHCND, months, percentiles, njitter, varname, stationCountsforQR, ...
			Tsubset, timeSubset, nboot, cacheDir);
		timeStamp = datestr(now);
		save(savename, 'QRpartial', 'timeStamp')
	end

	QR = [];

else
	% collect the results from the partial analyses

	if length(percentiles) == 3 % only do cross correlation if considering a small set of percentiles
		order = perms(1:length(percentiles));
		order = unique(order(:,1:2),'rows');
		order(order(:,1)>order(:,2),:) = [];
		OPTS.order = order;
		rho = NaN(nstations, size(OPTS.order,1));
	else
		rho = [];
	end
	
	QR.percentiles = percentiles;
	OPTS.varname = varname;
	OPTS.months = months;
	OPTS.njitter = njitter;
	OPTS.nboot = nboot;
	QR.opts = OPTS;

	beta = NaN(nstations, length(percentiles), 2);
	betaSD = beta;
	bootstrapSD = NaN(nstations, length(percentiles));
	bootstrapSD_OLS = NaN(nstations, 1);
	fragSig = NaN(nstations, length(percentiles));
	dbeta = NaN(nstations, 3);
	dbootstrapSD = NaN(nstations, 3);
	dfracSig = NaN(nstations, 3);

	savename0 = getHash([cacheDir '/QR'], varname, months, percentiles, njitter, ...
		nIntervals, Tsubset0, nboot);

	if exist(savename0, 'file')
		load(savename0);
	else

		for ct = 1:nIntervals
			stationCountsforQR = [intervalsQR(ct): intervalsQR(ct + 1)];
			timeInd = find(month(GHCND.time) >= min(months) & month(GHCND.time) <= max(months));
			Tsubset = GHCND.(varname)(timeInd, ismember(1:nstations, stationCountsforQR));
			timeSubset = GHCND.time(timeInd);
			savename = getHash([cacheDir '/QRstationSet'], varname, months, percentiles, njitter, Tsubset, ...
				stationCountsforQR, nboot);
 
			if ~exist(savename,'file');
				disp(['File does not exist for interval ' num2str(ct) '']);
			else
				load(savename);
				beta(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :, :) = ...
					QRpartial.beta;
				betaSD(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :, :) = ...
					QRpartial.betaSD;
				bootstrapSD(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
					QRpartial.bootstrapSD;
				bootstrapSD_OLS(ismember(1:nstations, QRpartial.opts.stationCountsforQR)) = ...
					QRpartial.bootstrapSD_OLS;
				dbeta(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
					QRpartial.dbeta;
				dbootstrapSD(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
					QRpartial.dbootstrapSD;
				dfracSig(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
					QRpartial.dfracSig;

				if length(percentiles) == 3
					rho(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
						QRpartial.rho;
				end

				fracSig(ismember(1:nstations, QRpartial.opts.stationCountsforQR), :) = ...
				 	QRpartial.fracSig;
			end

		end

		QR.beta = beta;
		QR.betaSD = betaSD;
		QR.bootstrapSD = bootstrapSD;
		QR.bootstrapSD_OLS = bootstrapSD_OLS;
		QR.rho = rho;
		QR.fracSig = fracSig;
		QR.dbeta = dbeta;
		QR.dbootstrapSD = dbootstrapSD;
		QR.dfracSig = dfracSig;

		timeStamp = datestr(now);
		save(savename0, 'QR', 'timeStamp');
	end

end