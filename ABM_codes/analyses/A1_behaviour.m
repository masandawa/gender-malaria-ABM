function A1_behaviour(X, m)
%A1_BEHAVIOUR  Effect of the female-male night-time activity contrast delta_A
%   on exposure and infection, and concentration of clinical and transmission
%   events (reference configuration, all other mechanisms on).
%   delta_A < 0: men more often outdoors. -0.15 reference (~16 min per night),
%   -1.0 occupational contrast (~1.8 h per night).
dGrid = [-1.0 -0.6 -0.3 -0.15 0 0.15 0.30];
T = struct([]); Cc = struct([]);
for r = 1:X.A1reps
    P = base_params(X, m); [P, C] = scenario(P, 'REF');
    set_seed(5000 + r); S0 = abm_init(P);
    S0 = abm_run(P, S0, X.burn, C);
    for id = 1:numel(dGrid)
        Pd = P; Pd.deltaA = dGrid(id);
        S = S0; [S.Pout, S.Bed] = activity_rows(Pd, S, (1:P.N)');
        set_seed(7000 + r);                                  % common random numbers
        [~, O] = abm_run(Pd, S, X.eval, C);
        R = summarize_run(Pd, O);
        T = addrow(T, 'rep', r, 'delta', dGrid(id), R);
        if dGrid(id) == -0.15
            Cc = addrow(Cc, 'rep', r, concentration(O));
        end
    end
    fprintf('[A1] replicate %d/%d\n', r, X.A1reps);
end
struct2csv('results/A1_behaviour.csv', T);
struct2csv('results/A1_concentration.csv', Cc);
end
