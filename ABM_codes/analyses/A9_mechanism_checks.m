function A9_mechanism_checks(X, m)
%A9_MECHANISM_CHECKS  Experiment 6: checks of the two interpretations on which
%   the Discussion relies.
%
%   Part 1 - clearance pathway and clinical immunity structure.
%     For three immunity structures (P.immType = 0 age-dependent reference;
%     1 infection-count immunity; 2 parasite-exposure immunity) the model is
%     re-calibrated to the same target and the equality, clearance-only (S) and
%     reference configurations are run from a common equality burn-in.
%     For people aged >= 5 years, by sex, it records the share of person-time
%     susceptible, the mean duration of non-clinical infection, new infections
%     and clinical episodes per person-year, and the mean clinical probability.
%
%   Part 2 - components of ITN protection (reference configuration).
%     ITNs are given to a random 70% of people at the start and kept for five
%     years (randomised receipt), and compared with a run without ITNs from the
%     same burn-in state. This separates (i) the exposure averted for a user,
%     (ii) the direct effect on holders (holders vs non-holders, randomised),
%     (iii) the indirect effect on non-holders (non-holders vs no ITNs) and
%     (iv) the population effect (all people, ITNs vs no ITNs).
%     The same runs give the simulated share of children <5 years and the
%     proportion of women aged 15-44 years who are pregnant (demography check).
if ~isfield(X, 'A9reps'), X.A9reps = 12; end
if ~isfield(X, 'A9calReps'), X.A9calReps = 4; end

% ------------------------------------------------------------- part 1
T = struct([]); M = struct([]);
cfgs = {[0 0 0 0], 'equality'; [0 0 0 1], 'clearance only'; [1 1 1 1], 'reference'};
for it = 0:2
    if it == 0
        mt = m;
    else
        P0 = base_params(X, m); P0.immType = it;
        mt = calibrate_m(X, P0, m, [0.8 1.0 1.25], X.A9calReps, 95000 + 100*it);
    end
    M = addrow(M, 'immType', it, 'm_star', mt);
    fprintf('[A9] immunity structure %d: m* = %.2f\n', it, mt);
    for r = 1:X.A9reps
        P = base_params(X, mt); P.immType = it;
        [Pb, Cb] = scenario(P, [0 0 0 0]);
        set_seed(91000 + r); S0 = abm_init(Pb);
        S0 = abm_run(Pb, S0, X.burn, Cb);
        for ic = 1:size(cfgs,1)
            [Pc, Cc] = scenario(P, cfgs{ic,1});
            S = S0; [S.Pout, S.Bed] = activity_rows(Pc, S, (1:P.N)');
            set_seed(92000 + r);
            [~, O] = abm_run(Pc, S, X.eval, Cc);
            R = summarize_run(Pc, O);
            T = addrow(T, 'immType', it, 'config', cfgs{ic,2}, 'rep', r, ...
                'prev', R.prev, 'inc', R.inc, 'incU5', R.incU5, 'PRfm', R.PRfm, 'IRRfm', R.IRRfm, ...
                pathway(O));
        end
    end
end
struct2csv('results/A9_immunity_runs.csv', T);
struct2csv('results/A9_immunity_calibration.csv', M);

% ------------------------------------------------------------- part 2
Q = struct([]);
for r = 1:X.A9reps
    P = base_params(X, m); [P, Cb] = scenario(P, 'REF');
    set_seed(93000 + r); S0 = abm_init(P);
    S0 = abm_run(P, S0, X.burn, Cb);
    Cn = struct('alloc', 'U', 'cov', 0.70, 'realloc', Inf);   % randomised, kept
    set_seed(94000 + r); [~, O1] = abm_run(P, S0, X.eval, Cn);
    C0 = struct('alloc', 'none', 'realloc', Inf);
    set_seed(94000 + r); [~, O0] = abm_run(P, S0, X.eval, C0);
    R1 = summarize_run(P, O1); R0 = summarize_run(P, O0);
    % per-person infections among people present throughout both runs
    stay = O1.born == O1.born0 & O0.born == O0.born0;
    z = O1.z0 == 1;
    yrs = X.eval/365;
    nH1 = mean(O1.Ninf(stay & z))/yrs;  nN1 = mean(O1.Ninf(stay & ~z))/yrs;
    nH0 = mean(O0.Ninf(stay & z))/yrs;  nN0 = mean(O0.Ninf(stay & ~z))/yrs;
    EX = O1.EX; eActU = sum(EX(:,2) - EX(:,4)); eEffU = sum(EX(:,9));
    Q = addrow(Q, 'rep', r, ...
        'avertedUser', 1 - eEffU/eActU, ...                 % share of a user's exposure averted
        'infHold', nH1, 'infNonHold', nN1, 'infNoITN_hold', nH0, 'infNoITN_non', nN0, ...
        'RRdirect', nH1/nN1, 'RRindirect', nN1/nN0, 'RRtotal', nH1/nH0, ...
        'prevITN', R1.prev, 'prevNoITN', R0.prev, 'incITN', R1.inc, 'incNoITN', R0.inc, ...
        'u5share', sum(O1.G(1:2,1))/sum(O1.G(:,1)), 'pregShare', O1.SP(1,1)/O1.SP(4,1));
    fprintf('[A9] ITN components replicate %d/%d\n', r, X.A9reps);
end
struct2csv('results/A9_itn_components.csv', Q);
end

% =========================================================================
function W = pathway(O)
% adults >= 5 years, columns F and M
ST = O.ST; G = O.G; py = G(3:4,1)/365;
W.suscF  = ST(1,1)/sum(ST(1,:));  W.suscM  = ST(2,1)/sum(ST(2,:));
W.infTF  = sum(ST(1,3:5))/sum(ST(1,:)); W.infTM = sum(ST(2,3:5))/sum(ST(2,:));
W.durIAF = ST(1,3)/max(O.IAX(1),1); W.durIAM = ST(2,3)/max(O.IAX(2),1);
W.infPYF = G(3,6)/py(1); W.infPYM = G(4,6)/py(2);
W.clinPYF = G(3,4)/py(1); W.clinPYM = G(4,4)/py(2);
W.pCF = O.PC(1,1)/max(O.PC(1,2),1); W.pCM = O.PC(2,1)/max(O.PC(2,2),1);
end
