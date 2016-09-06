function stationRecord = sampleHourly(USAF, WBAN, yearStart, yearEnd, timeOfMeasurement, csvDir, cacheDir)

% First, check to see if specified station has all data available
nyrs = length(yearStart:yearEnd);
files = dir([csvDir '/' (USAF) '-' (WBAN) '*.csv']);
types = '%s %s %d %d %d %d %d %s %f %f %f %f %s %f %s %f %s';
tempIndex = 12;
tempFlagIndex = 13;
  
if length(files) ~= nyrs
	disp(['Not all files are present for ' num2str(USAF) '-' num2str(WBAN) ''])
	stationRecord = [];
else  
	savename = getHash([cacheDir '/hourly2TxTn'], USAF, WBAN, yearStart, yearEnd, timeOfMeasurement, csvDir);
	if exist(savename, 'file')
		disp(['Loading file for ' num2str(USAF) '-' num2str(WBAN) ''])
		load(savename)
	else
		disp(['Calculating TMAX and TMIN for ' num2str(USAF) '-' num2str(WBAN) ''])
		% append values from each year
		TEMPERATURE = [];
		TIME = [];
		for ct = 1:nyrs
			fid = fopen([csvDir '/' files(ct).name]);

			headers = textscan(fid,'%q',17,'delimiter',',');
			vals = textscan(fid, types,'delimiter',',');
			fclose(fid);

			% turn 999.9 and bad flags to NaN
			% data that has passed all QC checks has a flag of 1, 5
			T = vals{tempIndex};
			T(T > 900) = NaN;
			flags = vals{tempFlagIndex};
			validValues = ismember(flags, '1') | ismember(flags, '5') | ismember(flags, '"1"') | ismember(flags, '"5"');
			T(~validValues) = NaN;
			time = datenum(double(vals{3}), double(vals{4}), double(vals{5}),...
				double(vals{6}), double(vals{7}), zeros(size(vals{7})));

			TEMPERATURE = [TEMPERATURE;T];
			TIME = [TIME;time];

		end

		% remove duplicate values for times
		% keep one if they have the same measurement, switch to NaN if they have different measurements
		d = diff(TIME);
		loc = find(d == 0);
		for ii = 1:length(loc)
			if TEMPERATURE(loc(ii)) ~= TEMPERATURE(loc(ii) + 1)
				TEMPERATURE(loc(ii)) = NaN;
			end
			TIME(loc(ii) + 1) = -999;
		end 
		pl = TIME == -999;
		TEMPERATURE(pl) = [];
		TIME(pl) = [];
 
		% for each day and time of measurement, calculate TMAX and TMIN for previous 24 hours
		dailyTime = datenum(yearStart,1,1):datenum(yearEnd,12,31);
		threshold = 1.5/24; % use value if measurement is within 1.5 hours of desired time (so 3 hrly data can be used)
		Tx = NaN(length(dailyTime), length(timeOfMeasurement));
		Tn = NaN(length(dailyTime), length(timeOfMeasurement));

		for jj = 1:length(dailyTime)
		
			for kk = 1:length(timeOfMeasurement)
		
				[val, indEnd] = min(abs(TIME - (dailyTime(jj) + timeOfMeasurement(kk)/24)));
				[val2, indStart] = min(abs(TIME - (dailyTime(jj) + timeOfMeasurement(kk)/24 - 1)));
				if val < threshold && val2 < threshold
					prev24 = TEMPERATURE(indStart:indEnd);
					Tx(jj,kk) = max(prev24);
					Tn(jj,kk) = min(prev24);
				end
				
			end
		end


		lat = vals{9}(1);
		lon = vals{10}(1);
		elev = vals{11}(1);

		stationRecord.Tx = Tx;
		stationRecord.Tn = Tn;
		stationRecord.loc = [lon lat elev];
		stationRecord.time = dailyTime;

		OPTS.USAF = USAF;
		OPTS.WBAN = WBAN;
		OPTS.timeOfMeasurement = timeOfMeasurement;
		OPTS.yearStart = yearStart;
		OPTS.yearEnd = yearEnd;
		stationRecord.opts = OPTS;

		timestamp = datestr(now);
		save(savename,'stationRecord','timestamp');

	end

end