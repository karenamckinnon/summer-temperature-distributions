function [qlevs0,ind_hmm,ind_imputed] = estimate_quantization_v2(t0,T0,varargin)

if nargin>2
    do_plot = varargin{1};
else
    do_plot = 0;
end
 
T_unique = unique(T0(~isnan(T0)));
T_unique_counts = sum( bsxfun(@eq, T0, T_unique') )';

% Evaluate possible values
S_C = [1 5 10];
S_f = [50/9];
S = [S_C S_f];
all_possible = (-700:700)'; % Almost-certainly-exhaustive range, -70C to 70C
P = NaN(length(all_possible),length(S));
for j = 1:length(S)
   if j<=length(S_C)
       possible = (all_possible(1):S(j):all_possible(end))';
   else
       possible = round(S(j)*(round((all_possible(1)/S(j) + 32):round(all_possible(end)/S(j) + 32))-32))';
   end
   P(:,j) = hist(possible,all_possible);
   P(:,j) = P(:,j) + 0*min(P(P(:,j)>0,j))/1000; % Tiny perturbation to avoid failing completely on transcription errors.
                                              % (This is hard coded, but
                                              % feel free to remove if you
                                              % feel the data are perfect
                                              % and that your list of
                                              % states is exhaustive).
   P(:,j) = P(:,j)./sum(P(:,j));
end
ind_hmm = ~isnan(T0) & T0<700 & T0>-400;
ind_imputed = ~ind_hmm;
Ts = T0(ind_hmm);

O = unique(Ts); % Search only over the observed subset -- faster without loss of generality.
% Subset P
seq = Ts;
Ps = NaN(length(O),length(S));
for j = 1:length(O)
    seq(Ts==O(j)) = j; % Sequence is just a timeseries of observation index labels.
    ind = find(all_possible == O(j));
    Ps(j,:) = P(ind,:);% Collapse the MC estimate of the priors for use in the emission matrix.
end 

alpha = 0.999; % Hard-coded, very small probability of changing states.
trans = eye(length(S))*alpha + ... % Transition matrix (assume equal distribution in the non-persistent transitions).
        triu((1-alpha)/(length(S)-1)*ones(length(S)),1) + ...
        tril((1-alpha)/(length(S)-1)*ones(length(S)),-1);
    
emis = Ps'; % Emission probability matrix (numstates x numobservables)
estimatedStates = hmmviterbi(seq,trans,emis); % Solve the hidden Markov model.
qlevs = S(estimatedStates); % Actual numerical quantization estimates, pointwise.
% Fill in non-estimated states by nearest-neighbor (these are likely NaN
% values in T0, but the estimates are provided anyways).
if sum(ind_imputed) > 0
    idx = knnsearch(find(ind_hmm),find(ind_imputed));
    qlevs0 = NaN(size(T0));
    qlevs0(ind_hmm) = qlevs;
    qlevs0(ind_imputed) = qlevs(idx);
else
    qlevs0 = qlevs;
end

if do_plot
    figure(2)
    clf
    h1=subplot(3,1,1);
    plot(t0(ind_hmm),T0(ind_hmm),'b.')
    hold on
    plot(t0(ind_imputed),T0(ind_imputed),'rx')
    hold off
    ylabel('T_{max} (C)')
    datetick
    h2=subplot(3,1,2);
    plot(t0(ind_hmm),T0(ind_hmm),'b.')
    hold on
    plot(t0(ind_imputed),T0(ind_imputed),'rx')
    hold off
    ylabel('\Delta T (C)')
    datetick
    h3=subplot(3,1,3);
    plot(t0(ind_hmm),qlevs0(ind_hmm),'bo')
    hold on
    plot(t0(ind_imputed),qlevs0(ind_imputed),'rx')
    hold off
    grid on
    ylabel('Estimated Quantization')
    datetick
    linkaxes([h1 h2 h3],'x')
end

return