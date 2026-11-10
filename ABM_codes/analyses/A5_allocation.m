function A5_allocation(X, m)
%A5_ALLOCATION  Individual-level ITN assignment experiment with a fixed number
%   of nets (70% of the population), re-selected annually:
%     U uniform; P policy-aligned (pregnant women, then children <5, then random);
%     D age-sex groups ranked by clinical incidence in the previous year;
%     E ideal-information expected exposure (benchmark, not observable);
%     A observable activity proxy (age, pregnancy, night-time activity, no e0);
%     H infections in the previous year; T attributed transmission in the
%     previous year (benchmark, not observable).
%   This is not a simulation of household-based national ITN distribution.
rules = {'U', 'P', 'D', 'E', 'A', 'H', 'T'};
T = struct([]);
for r = 1:X.A5reps
    P = base_params(X, m); [P, Cb] = scenario(P, 'REF');
    set_seed(20000 + r); S0 = abm_init(P);
    S0 = abm_run(P, S0, X.burn, Cb);
    Ck = Cb; Ck.alloc = 'keep';
    [S0, Ob] = abm_run(P, S0, X.hist, Ck);            % history window for D, H and T
    J = Ob.G(:,4)./Ob.G(:,1); [~, Dorder] = sort(J, 'descend');
    for is = 1:numel(rules)
        C = struct('alloc', rules{is}, 'cov', 0.70, 'realloc', 365, 'Dorder', Dorder');
        set_seed(30000 + r);
        [~, O] = abm_run(P, S0, X.eval, C);
        T = addrow(T, 'rep', r, 'rule', rules{is}, summarize_run(P, O));
    end
    fprintf('[A5] replicate %d/%d\n', r, X.A5reps);
end
struct2csv('results/A5_allocation.csv', T);
end
