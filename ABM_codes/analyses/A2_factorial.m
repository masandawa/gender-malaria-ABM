function A2_factorial(X, m)
%A2_FACTORIAL  Full 2^4 factorial experiment over the four mechanism blocks
%   B (gendered behaviour), A (programme access), G (pregnancy biology) and
%   S (sex-specific clearance); see model/scenario.m. Each replicate starts
%   from a common equality burn-in and all 16 configurations share random
%   numbers. Main effects, interactions and Shapley contributions are
%   computed in plotting/make_tables.m.
T = struct([]);
cfg = dec2bin(0:15) == '1';                         % rows: [B A G S]
for r = 1:X.A2reps
    P = base_params(X, m);
    [Pb, Cb] = scenario(P, [0 0 0 0]);
    set_seed(9000 + r); S0 = abm_init(Pb);
    S0 = abm_run(Pb, S0, X.burn, Cb);
    for ic = 1:16
        [Pc, Cc] = scenario(P, cfg(ic,:));
        S = S0; [S.Pout, S.Bed] = activity_rows(Pc, S, (1:P.N)');
        set_seed(11000 + r);
        [~, O] = abm_run(Pc, S, X.eval, Cc);
        R = summarize_run(Pc, O);
        T = addrow(T, 'rep', r, 'B', cfg(ic,1), 'A', cfg(ic,2), 'G', cfg(ic,3), 'S', cfg(ic,4), R);
    end
    fprintf('[A2] replicate %d/%d\n', r, X.A2reps);
end
struct2csv('results/A2_factorial.csv', T);
end
