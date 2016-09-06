function GHCND = getGHCND(varname, ghcnd_matpath, yearStart, yearEnd, months, latRange, lonRange, cacheDir);
 
savename = getHash([cacheDir '/ghcndData'], varname, ghcnd_matpath, yearStart, yearEnd, ...
	months, latRange, lonRange);

if exist(savename)
	load(savename)
else

	% USER NEEDS TO HAVE DATABASE OF GHCND DATA -- CONTACT AUTHORS FOR MORE INFORMATION
	% A sample metadata.mat and station matfile are included in the repo
	% load metadata
	load([ghcnd_matpath '/metadata.mat']) 
	fn = char(['metadata.' varname '_nobs']);
	fn1 = char(['metadata.' varname '_firstdate']);
	fn2 = char(['metadata.' varname '_lastdate']);
	fnlat = char('metadata.lat');
	fnlon = char('metadata.lon');
	date1 = datevec(eval(fn1));
	date2 = datevec(eval(fn2));

	% pull out desired region of interest
	index =  eval(fnlat) > min(latRange) & eval(fnlat) < max(latRange) & ...
		eval(fnlon) > min(lonRange) & eval(fnlon) < max(lonRange);

	files_with_data = find(eval(fn) > 0 & date1(:,1)' <= yearStart & date2(:,1)' >= yearEnd & index);

	numyr = length(yearStart:yearEnd);
	% initialize all matrices
	nf = length(files_with_data);

	TIME = datenum(yearStart, months(1), 1):datenum(yearEnd, months(end), eomday(yearEnd, months(end)));
	GHCND.loc = NaN(nf, 3);
	GHCND.(varname) = NaN(length(TIME), nf);
	GHCND.id = NaN(nf, 11);

	for ct = 1:nf  

		fn = [ghcnd_matpath '/' (metadata.station_codes{files_with_data(ct)})];
		SR = load(fn);
		stationRecord = SR.stationRecord;
		disp(['Station: ' stationRecord.meta.id ', ' num2str(ct) ' out of ' num2str(nf) ''])

		% Get temperature and date vector. flag = 1 if fits date criteria
		[T d flag] = getGHCNDfile(yearStart,yearEnd,varname,stationRecord);

		% Make sure that record spans months of interest (MAY WANT TO CHANGE LATER)
		if year(d(1)) == yearStart && month(d(1)) == months(1) && day(d(1)) ~= 1,
			flag = 0;
		end
		if year(d(end)) == yearEnd && month(d(end)) == months(2) && day(d(end)) ~= eomday(yearEnd, months(2)),
			flag = 0;
		end

		if flag
			% pull out desired years
			% Add NaNs if part of the year is missing
			GHCND.loc(ct,:) = [stationRecord.meta.lon, stationRecord.meta.lat, stationRecord.meta.elev];
			indexUse = year(d) >= yearStart & year(d) <= yearEnd & month(d) >= months(1) & month(d) <= months(end);
			Tdummy = NaN(size(TIME));
			timeUse = d(indexUse);
			Tdummy(ismember(TIME, timeUse)) = T(indexUse);
			GHCND.(varname)(:, ct) = Tdummy;
			GHCND.id(ct, :) = stationRecord.meta.id;
		end

	end

	GHCND.time = TIME;
	GHCND.id = char(GHCND.id);

	timestamp = datestr(now);
	inputs.varname = varname;
	inputs.ghcnd_matpath = ghcnd_matpath;
	inputs.yearStart = yearStart;
	inputs.yearEnd = yearEnd;
	inputs.months = months;
	inputs.lonRange = lonRange;
	inputs.latRange = latRange;

	save(savename, 'GHCND', 'timestamp', 'inputs')

end

return