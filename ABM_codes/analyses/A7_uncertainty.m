function A7_uncertainty(X, m, draws)
%A7_UNCERTAINTY  Joint parameter uncertainty (Latin hypercube sampling).
%   Nine uncertain parameters are sampled jointly from uniform ranges. For each
%   draw the model is run once for: ITN selection rules P, U, D, E and A
%   (reference configuration), and the single mechanisms "sex-specific
%   clearance only" and "occupational activity contrast only". The analysis
%   reports how often each principal ranking or mechanism comparison holds.
%   Stochastic and parameter uncertainty are combined (one realisation/draw).
%   For every draw the mosquito-to-human ratio m is re-calibrated to the same
%   under-five prevalence target (calibrate_m), so that the sampled parameters
%   change the mechanisms but not baseline transmission intensity. The
%   re-calibrated value is stored as m_cal.
pr = {'dA', 150, 250; 'etaN', 0.50, 0.90; 'hrClear', 1.20, 2.75; ...
      'kappaP', 1.5, 2.5; 'rrClinP', 1.0, 2.0; 'sigAct', 0.25, 0.75; ...
      'Rout', 0.70, 1.43; 'qA', 0.05, 0.20};
n = X.A7draws; k = size(pr,1);
set_seed(60000);
U = zeros(n,k);
for j = 1:k, U(:,j) = (randperm(n)' - rand(n,1))/n; end   % Latin hypercube on (0,1)
if nargin < 3, draws = 1:n; end
T = struct([]);
for i = draws
    P = base_params(X, m);
    lab = struct();
    for j = 1:k
        v = pr{j,2} + U(i,j)*(pr{j,3} - pr{j,2}); P.(pr{j,1}) = v; lab.(pr{j,1}) = v;
    end
    P.m = calibrate_m(X, P, m, [0.7 1.0 1.4], X.A7calReps, 66000 + 100*i);
    lab.m_cal = P.m;
    row = lab;
    % ---- ITN selection rules, reference configuration
    [Pr, Cr] = scenario(P, 'REF');
    set_seed(61000 + i); S0 = abm_init(Pr);
    S0 = abm_run(Pr, S0, X.burn, Cr);
    Ck = Cr; Ck.alloc = 'keep';
    [S0, Ob] = abm_run(Pr, S0, X.hist, Ck);
    J = Ob.G(:,4)./Ob.G(:,1); [~, Dorder] = sort(J, 'descend');
    for s = {'P', 'U', 'D', 'E', 'A'}
        C = struct('alloc', s{1}, 'cov', 0.70, 'realloc', 365, 'Dorder', Dorder');
        set_seed(62000 + i);
        [~, O] = abm_run(Pr, S0, X.eval, C); R = summarize_run(Pr, O);
        row.(['prev_' s{1}]) = R.prev; row.(['inc_' s{1}]) = R.inc;
        row.(['preg_' s{1}]) = R.rdtPreg; row.(['incU5_' s{1}]) = R.incU5;
    end
    % ---- single mechanisms from equality
    [Pb, Cb] = scenario(P, [0 0 0 0]);
    set_seed(63000 + i); S0 = abm_init(Pb);
    S0 = abm_run(Pb, S0, X.burn, Cb);
    [Ps, Cs] = scenario(P, [0 0 0 1]);
    set_seed(64000 + i); [~, O] = abm_run(Ps, S0, X.eval, Cs); R = summarize_run(Ps, O);
    row.PR_clear = R.PRfm; row.dPrev_clear = R.dPrev;
    Po = Pb; Po.deltaA = -1.0;
    S = S0; [S.Pout, S.Bed] = activity_rows(Po, S, (1:P.N)');
    set_seed(64000 + i); [~, O] = abm_run(Po, S, X.eval, Cb); R = summarize_run(Po, O);
    row.PR_occ = R.PRfm; row.dPrev_occ = R.dPrev;
    if isempty(T), T = row; else, T(end+1) = orderfields(row, T(1)); end %#ok<AGROW>
    fprintf('[A7] draw %d/%d\n', i, n);
    if ~isequal(draws, 1:n)              % partial runs save after every draw
        struct2csv(sprintf('results/A7_uncertainty_%d_%d.csv', draws(1), draws(end)), T);
    end
end
if isequal(draws, 1:n)
    struct2csv('results/A7_uncertainty.csv', T);
else                                   % partial run for parallel use; rows are identical
    struct2csv(sprintf('results/A7_uncertainty_%d_%d.csv', draws(1), draws(end)), T);
end
end
