function [ISD] = sampleHourlyWrapper(allUSAF, allWBAN, timesToSample, runNumber, ...
    nIntervals, opts)

% Get station count
nstations = size(allUSAF, 1); 
intervals = round(linspace(1, nstations, nIntervals + 1));
if intervals(end) ~= nstations; intervals(end) = nstations; end

if runNumber > 0

    stationCounts = [(intervals(runNumber)): intervals(runNumber + 1)];
 
    savename = getHash([opts.cacheDir '/hourlySamplesSet'], runNumber, nIntervals, allUSAF, ...
        allWBAN, opts.yearStart, opts.yearEnd, timesToSample);

    if exist(savename)
        load(savename)
        disp(['Loading file set ' num2str(runNumber) ''])
    else
        counter = 1;
        for ct = stationCounts
        	disp(['Starting station ' num2str(ct) ' of ' num2str(stationCounts(end)) ''])
            tic
            stationRecord = sampleHourly(allUSAF(ct, :), allWBAN(ct, :), opts.yearStart, opts.yearEnd, ...
                timesToSample, opts.csvDir, opts.cacheDir);

            if ~isempty(stationRecord)
                ISD.loc(counter, :) = stationRecord.loc;
                ISD.TMAX(:, :, counter) = stationRecord.Tx;
                ISD.TMIN(:, :, counter) = stationRecord.Tn;
                ISD.time = stationRecord.time;
                ISD.USAF(counter, :) = stationRecord.opts.USAF;
                ISD.WBAN(counter, :) = stationRecord.opts.WBAN;

                ISD.opts = opts;
        		ISD.opts.timeOfMeasurement = stationRecord.opts.timeOfMeasurement;

                counter = counter + 1;
            end
            toc

        end


        ISDpartial = ISD;
        clear ISD

        timestamp = datestr(now);
        save(savename, 'ISDpartial', 'timestamp' ,'opts')
    end

    ISD = [];
 

else

    savename0 = getHash([opts.cacheDir '/hourly2TxTn'], nIntervals, allUSAF, ...
        allWBAN, opts.yearStart, opts.yearEnd, timesToSample);

    if exist(savename0, 'file')
        load(savename0);
    else

        ISD.loc = [];
        ISD.TMAX = [];
        ISD.TMIN = [];
        ISD.USAF = [];
        ISD.WBAN = [];


        for ct = 1:nIntervals
            stationCounts = [intervals(ct): intervals(ct + 1)];
            savename = getHash([opts.cacheDir '/hourlySamplesSet'], ct, nIntervals, allUSAF, ...
                allWBAN, opts.yearStart, opts.yearEnd, timesToSample);

            if ~exist(savename,'file');
                disp(['File does not exist for interval ' num2str(ct) '']);
            else
                load(savename)

                ISD.loc = cat(1, ISD.loc, ISDpartial.loc);
                ISD.TMAX = cat(3, ISD.TMAX, ISDpartial.TMAX);
                ISD.TMIN = cat(3, ISD.TMIN, ISDpartial.TMIN);
                ISD.USAF = cat(1, ISD.USAF, ISDpartial.USAF);
                ISD.WBAN = cat(1, ISD.WBAN, ISDpartial.WBAN);
                ISD.time = ISDpartial.time;
                ISD.opts = ISDpartial.opts;
            end

            timeStamp = datestr(now);
            save(savename0, 'ISD', 'timeStamp');
        end
    end

end

