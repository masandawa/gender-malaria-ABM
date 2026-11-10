function [mt, b] = calibrate_m(X, P0, m0, fac, reps, seed0)
%CALIBRATE_M  Re-calibrate the mosquito-to-human ratio m to the under-five RDT
%   target (P0.calTarget) for parameter set P0 (reference configuration).
%   m is evaluated at m0*fac with 'reps' realisations each (burn-in X.burn,
%   one measurement year); m* inverts a log-log regression of prevalence on m,
%   as in A0_calibration.
mv = []; yv = [];
for ig = 1:numel(fac)
    for k = 1:reps
        P = P0; P.m = m0*fac(ig); [P, C] = scenario(P, 'REF');
        set_seed(seed0 + 10*ig + k); S = abm_init(P);
        S = abm_run(P, S, X.burn, C);
        [~, O] = abm_run(P, S, 365, C); R = summarize_run(P, O);
        mv(end+1) = P.m; yv(end+1) = R.u5rdt; %#ok<AGROW>
    end
end
ok = yv > 0;
b = [ones(sum(ok),1) log(mv(ok)')] \ log(yv(ok)');
mt = exp((log(P0.calTarget) - b(1))/b(2));
mt = min(max(mt, 0.3*m0), 3*m0);     % guard against extrapolation
end
