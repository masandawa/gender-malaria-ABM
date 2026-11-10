function A6_biting_ratio(X, m)
%A6_BITING_RATIO  Sensitivity to the outdoor:indoor biting-rate ratio R_out/in.
%   Two questions are answered separately:
%     'fixed' : mosquito density held at the reference m* (biting location and
%               transmission intensity both change);
%     'recal' : m re-calibrated to the same under-five prevalence target at each
%               ratio (only the location of biting changes).
%   For each ratio: behavioural female-male prevalence differences under
%   equality of all other mechanisms (delta_A = 0, -0.15, -1.0), and the
%   comparison of policy-aligned (P) and ideal-information exposure (E) ITN
%   selection in the reference configuration.
Rgrid = unique([1.00 X.Rgrid]);
T = struct([]); Cal = struct([]);
for iR = 1:numel(Rgrid)
    Rv = Rgrid(iR);
    modes = {'fixed', 'recal'}; if Rv == 1, modes = {'fixed'}; end
    for im = 1:numel(modes)
        mv = m;
        if strcmp(modes{im}, 'recal')
            mv = recalibrate(X, m, Rv);
            Cal = addrow(Cal, 'Rout', Rv, 'm_star', mv);
        end
        for r = 1:X.A6reps
            P = base_params(X, mv); P.Rout = Rv;
            % behavioural contrasts with all other mechanisms at equality
            [Pb, Cb] = scenario(P, [0 0 0 0]);
            set_seed(40000 + r); S0 = abm_init(Pb);
            S0 = abm_run(Pb, S0, X.burn, Cb);
            for d = [0 -0.15 -1.0]
                Pd = Pb; Pd.deltaA = d;
                S = S0; [S.Pout, S.Bed] = activity_rows(Pd, S, (1:P.N)');
                set_seed(42000 + r);
                [~, O] = abm_run(Pd, S, X.eval, Cb);
                T = addrow(T, 'Rout', Rv, 'mode', modes{im}, 'm', mv, 'rep', r, ...
                    'run', sprintf('delta%.2f', d), summarize_run(Pd, O));
            end
            % allocation comparison in the reference configuration
            [Pr, Cr] = scenario(P, 'REF');
            set_seed(44000 + r); S0 = abm_init(Pr);
            S0 = abm_run(Pr, S0, X.burn, Cr);
            for s = {'P', 'E'}
                C = struct('alloc', s{1}, 'cov', 0.70, 'realloc', 365);
                set_seed(46000 + r);
                [~, O] = abm_run(Pr, S0, X.eval, C);
                T = addrow(T, 'Rout', Rv, 'mode', modes{im}, 'm', mv, 'rep', r, ...
                    'run', ['rule' s{1}], summarize_run(Pr, O));
            end
        end
        fprintf('[A6] R_out/in = %.2f (%s) done\n', Rv, modes{im});
    end
end
struct2csv('results/A6_biting_ratio.csv', T);
if ~isempty(Cal), struct2csv('results/A6_recalibration.csv', Cal); end
end

function mv = recalibrate(X, m, Rv)
% re-calibrate m to the same under-five target at ratio Rv
P = base_params(X, m); S = abm_init(P);
om = 1 - P.rho*exp(-S.age/P.aage);
f1 = mean(S.e0.*om.*(S.Pout*P.w' + (1-S.Pout)*P.w'));
fR = mean(S.e0.*om.*(S.Pout*(P.w*Rv)' + (1-S.Pout)*P.w'));
grid = m*(f1/fR)*[0.8 1.0 1.25]; mm = []; yy = [];
for ig = 1:3
    for k = 1:max(2, round(X.calReps/2))
        P = base_params(X, grid(ig)); P.Rout = Rv; [P, C] = scenario(P, 'REF');
        set_seed(48000 + 10*ig + k); S = abm_init(P);
        S = abm_run(P, S, X.burn, C); [~, O] = abm_run(P, S, 365, C);
        R = summarize_run(P, O); mm(end+1) = grid(ig); yy(end+1) = R.u5rdt; %#ok<AGROW>
    end
end
ok = yy > 0; b = [ones(sum(ok),1) log(mm(ok)')] \ log(yy(ok)');
mv = exp((log(P.calTarget) - b(1))/b(2));
end
