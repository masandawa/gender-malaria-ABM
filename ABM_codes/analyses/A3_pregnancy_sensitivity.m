function A3_pregnancy_sensitivity(X, m)
%A3_PREGNANCY_SENSITIVITY  One-at-a-time sensitivity of pregnancy outcomes to
%   the three uncertain pregnancy parameters (reference configuration):
%   attractiveness kappa_P, clinical multiplier r_C^P, IPTp efficacy.
%   Outputs compare pregnant with non-pregnant women aged 15-44 years.
V = {'reference', '', NaN;
     'kappa_P = 1.0', 'kappaP', 1.0; 'kappa_P = 1.5', 'kappaP', 1.5; 'kappa_P = 2.5', 'kappaP', 2.5;
     'r_C = 1.0', 'rrClinP', 1.0;   'r_C = 2.0', 'rrClinP', 2.0;
     'IPTp = 0.2', 'effIPTp', 0.2;   'IPTp = 0.6', 'effIPTp', 0.6};
T = struct([]);
for r = 1:X.A3reps
    P = base_params(X, m); [P, C] = scenario(P, 'REF');
    set_seed(13000 + r); S0 = abm_init(P);
    S0 = abm_run(P, S0, X.burn, C);
    for iv = 1:size(V,1)
        Pv = P; if ~isempty(V{iv,2}), Pv.(V{iv,2}) = V{iv,3}; end
        set_seed(15000 + r);
        [~, O] = abm_run(Pv, S0, X.eval, C);
        T = addrow(T, 'rep', r, 'variant', V{iv,1}, summarize_run(Pv, O));
    end
    fprintf('[A3] replicate %d/%d\n', r, X.A3reps);
end
struct2csv('results/A3_pregnancy.csv', T);
end
