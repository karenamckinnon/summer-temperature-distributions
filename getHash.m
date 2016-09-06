function fname = getHash(baseFname, varargin);
% wrapper for DataHash to check if file already exists based on inputs

rng('default')
for ct = 1:(nargin - 1) 
	cache_input.(['input' num2str(ct) '']) = varargin{ct};
end

cacheName = DataHash(cache_input);

fname = [baseFname '_' cacheName '.mat'];
