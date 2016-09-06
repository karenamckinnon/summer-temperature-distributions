function GHCND = doQC(GHCND, months, varname)
% Records from individual weather stations are included only if they have at least
% 80% coverage during desired months for at least 80% of the years considered
% in the analysis.
 
% NOTE THAT MONTHS MUST BE IN THE SAME YEAR FOR THE CURRENT APPROACH, i.e will not work for DJF

timeInd = find(month(GHCND.time) >= min(months) & month(GHCND.time) <= max(months));
Tsubset = GHCND.(varname)(timeInd, :);
timeSubset = GHCND.time(timeInd);
  
% remove leap days for this calculation
if ismember(2, months),
	Tsubset(month(timeSubset) == 2 & day(timeSubset) == 29, :) = [];
	timeSubset(month(timeSubset) == 2 & day(timeSubset) == 29) = [];
end

nyrs = length(unique(year(GHCND.time(timeInd))));
ndays = size(Tsubset, 1)/nyrs;
nstations = size(Tsubset, 2);

temp_mat = reshape(Tsubset, [ndays nyrs nstations]);
% fraction missing per set of months (station x year)
frac_missing = squeeze(sum(isnan(temp_mat), 1))/ndays;
% 80th percentile across years of the fraction of missing days
frac_missing_cutoff = prctile(frac_missing, 80, 1);
% only include stations where the 80th percentile is less than 0.2
% i.e. at least 80% of years must have 80% of data
idx_include = find(frac_missing_cutoff < 0.2);
GHCND.(varname) = GHCND.(varname)(:, idx_include);
GHCND.loc = GHCND.loc(idx_include, :);
if isfield(GHCND, 'id')
	GHCND.id = GHCND.id(idx_include, :);
end

