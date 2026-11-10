function A4_itn_use(X, m)
%A4_ITN_USE  Realised ITN protection under habitual (reference) and
%   independent nightly use, and whether the comparison of policy-aligned (P)
%   and exposure-informed (E) ITN selection depends on this choice.
%   Protection is summarised by the number of new infections per person during
%   follow-up among ITN holders relative to non-holders, directly standardised
%   over quintiles of baseline expected exposure (rate ratio), and by five-year
%   cumulative infection risk standardised in the same way.
T = struct([]);
for r = 1:X.A4reps
    for nightly = [0 1]
        P = base_params(X, m); P.useNightly = nightly; [P, Cb] = scenario(P, 'REF');
        set_seed(17000 + r); S0 = abm_init(P);
        S0 = abm_run(P, S0, X.burn, Cb);
        set_seed(19000 + r);
        [~, O] = abm_run(P, S0, X.eval, Cb);
        T = addrow(T, 'rep', r, 'nightly', nightly, 'rule', 'policyProb', summarize_run(P, O), rate_ratio(O));
        for s = {'P', 'E'}
            C = struct('alloc', s{1}, 'cov', 0.70, 'realloc', 365);
            set_seed(19000 + r);
            [~, O] = abm_run(P, S0, X.eval, C);
            T = addrow(T, 'rep', r, 'nightly', nightly, 'rule', s{1}, summarize_run(P, O), rate_ratio(O));
        end
    end
    fprintf('[A4] replicate %d/%d\n', r, X.A4reps);
end
struct2csv('results/A4_itn_use.csv', T);
end

function Q = rate_ratio(O)
% exposure-standardised infection rate among ITN holders vs non-holders
stay = O.born == O.born0; e = O.eExp0(stay); n = O.Ninf(stay); z = O.z0(stay) == 1;
[~, o] = sort(e); q = zeros(numel(e),1); q(o) = min(5, 1 + floor(5*(0:numel(e)-1)'/numel(e)));
rr = NaN(5,1); rn = NaN(5,1);
for k = 1:5, rr(k) = mean(n(q == k & z)); rn(k) = mean(n(q == k & ~z)); end
ok = isfinite(rr) & isfinite(rn);
Q.infRec = mean(rr(ok)); Q.infNon = mean(rn(ok)); Q.IRRinf = Q.infRec/max(Q.infNon, eps);
end
