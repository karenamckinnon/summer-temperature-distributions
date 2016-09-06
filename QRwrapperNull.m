function QRnull = QRwrapperNull(GHCND, varname, nIntervals, months, runNumber, percentiles, ...
	cacheDir, nboot);
 
% Get station count
nstations = size(GHCND.(varname), 2);
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

	savename = getHash([cacheDir '/QRnullstationSet'], months, percentiles, varname, ...
		nboot, timeSubset, Tsubset);

	if exist(savename, 'file')
		load(savename)
	else
		QRpartial = doQRnull(GHCND, months, percentiles, varname, nboot, timeSubset, ...
			Tsubset, cacheDir);
		timeStamp = datestr(now);
		save(savename, 'QRpartial', 'timeStamp')
	end

	QRnull = []; % don't output anything -- just doing the calculations

else
	% collect different station counts

	QRnull.percentiles = percentiles;
	OPTS.varname = varname;
	OPTS.months = months;
	OPTS.nboot = nboot;
	QRnull.opts = OPTS;
 
	betaNull = NaN(nstations, length(percentiles), nboot);

	savename0 = getHash([cacheDir '/QRnull'], varname, months, percentiles, ...
		nIntervals, Tsubset0, nboot);

	if exist(savename0, 'file')
		load(savename0);
	else

		for ct = 1:nIntervals
			stationCountsforQR = [intervalsQR(ct): intervalsQR(ct + 1)];
			timeInd = find(month(GHCND.time) >= min(months) & month(GHCND.time) <= max(months));
			Tsubset = GHCND.(varname)(timeInd, ismember(1:nstations, stationCountsforQR));
			timeSubset = GHCND.time(timeInd);


			savename = getHash([cacheDir '/QRnullstationSet'], months, percentiles, varname, ...
				nboot, timeSubset, Tsubset);

			if ~exist(savename,'file');
				disp(['File does not exist for interval ' num2str(ct) '']);
			else
				load(savename);
				betaNull(ismember(1:nstations, stationCountsforQR), :, :) = ...
					QRpartial.betaNull;

			end

		end

		QRnull.betaNull = betaNull;


		timeStamp = datestr(now);
		save(savename0, 'QRnull', 'timeStamp');
	end

end