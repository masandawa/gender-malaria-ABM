function make_tables()
%MAKE_TABLES  Summarise results/*.csv into the numbers reported in the paper.
%   Writes results/summary_*.csv (one file per table or figure) and
%   results/summary.txt. Means are over realisations; intervals are 95%
%   Monte Carlo intervals (t-based); paired comparisons use common random numbers.
fid = fopen('results/summary.txt', 'w');
cleanup = onCleanup(@() fclose(fid));

% ------------------------------------------------------------ A0 calibration
if exist('results/A0_calibration_summary.csv', 'file')
    S = csv2struct('results/A0_calibration_summary.csv');
    fprintf(fid, 'A0 calibration: m* = %.2f (bootstrap %.2f-%.2f; acceptance %.2f-%.2f; slope %.2f)\n', ...
        S.m_star, S.boot_lo, S.boot_hi, S.acc_lo, S.acc_hi, S.slope);
    Pf = csv2struct('results/A0_profile.csv');
    for i = 1:numel(Pf)
        fprintf(fid, '   profile %-14s m* %.2f  incidence %.0f  prevalence %.3f\n', Pf(i).set, Pf(i).m_star, Pf(i).inc, Pf(i).prev);
    end
end

% ------------------------------------------------------------ empirical evaluation (A2 reference)
if exist('results/A2_factorial.csv', 'file')
    F = csv2struct('results/A2_factorial.csv');
    ref = F([F.B] == 1 & [F.A] == 1 & [F.G] == 1 & [F.S] == 1);
    ev = {'u5rdt', 'RDT prevalence 6-59 months'; 'piI', 'Indoor share of exposure, non-users'; ...
          'piIF', 'Indoor share, females'; 'piIM', 'Indoor share, males'; 'use', 'ITN use, all'; ...
          'effCovPreg', 'ITN use, pregnant women'; 'PRu5', 'F:M prevalence <5 y'; ...
          'PRfmRDT', 'F:M RDT prevalence, all ages'; 'PRfm', 'F:M prevalence (all infections)'; ...
          'PRadult', 'F:M prevalence 15-44 y'; 'IRRfm', 'F:M clinical incidence'; ...
          'PregU5', 'Pregnant : children 6-59 m (RDT)'; 'ORu5', 'OR malaria, <5 non-users vs users'; ...
          'rdt', 'RDT prevalence, all ages'; 'prev', 'Prevalence, all infections'; 'inc', 'Clinical incidence /1000 PY'};
    E = struct([]);
    fprintf(fid, '\nEmpirical evaluation (reference configuration)\n');
    for i = 1:size(ev,1)
        q = ci([ref.(ev{i,1})]);
        E = addrow(E, 'quantity', ev{i,1}, 'mean', q(1), 'lo', q(2), 'hi', q(3));
        fprintf(fid, '   %-38s %.3f (%.3f-%.3f)\n', ev{i,2}, q);
    end
    struct2csv('results/summary_evaluation.csv', E);
end

% ------------------------------------------------------------ A1 behaviour
if exist('results/A1_behaviour.csv', 'file')
    B = csv2struct('results/A1_behaviour.csv'); d = [B.delta]; reps = unique([B.rep]);
    dg = unique(d); A1 = struct([]);
    fprintf(fid, '\nA1 behaviour (reference configuration)\n');
    for x = dg
        k = abs(d - x) < 1e-9; k0 = abs(d) < 1e-9;
        beh = [B(k).dPrev] - [B(k0).dPrev];            % behavioural component, paired
        q = ci(beh); p = ci([B(k).PRfm]);
        A1 = addrow(A1, 'delta', x, 'dBeh', q(1), 'dBeh_lo', q(2), 'dBeh_hi', q(3), ...
            'PRfm', p(1), 'outF', mean([B(k).outF]), 'outM', mean([B(k).outM]), ...
            'piIF', mean([B(k).piIF]), 'piIM', mean([B(k).piIM]));
        fprintf(fid, '   delta %+5.2f behavioural F-M difference %+.3f pp (%+.3f,%+.3f); PR %.3f; outdoor F %.3f M %.3f\n', ...
            x, q, p(1), mean([B(k).outF]), mean([B(k).outM]));
    end
    struct2csv('results/summary_A1.csv', A1);
    Cc = csv2struct('results/A1_concentration.csv');
    f = {'pClin','pTrans','qB50','qB80','qT50','qT80','qBt50','qTt50','J','meanInf','inf0','inf1','inf2','inf3p'};
    Cs = struct();
    for i = 1:numel(f), q = ci([Cc.(f{i})]); Cs.(f{i}) = q(1); Cs.([f{i} '_lo']) = q(2); Cs.([f{i} '_hi']) = q(3);
        fprintf(fid, '   %-8s %.4f (%.4f-%.4f)\n', f{i}, q); end
    struct2csv('results/summary_A1_concentration.csv', Cs);
end

% ------------------------------------------------------------ A2 factorial and Shapley
if exist('results/A2_factorial.csv', 'file')
    fprintf(fid, '\nA2 factorial (16 configurations)\n');
    outs = {'logPR', 'logIRR', 'dEff', 'Gabs'};
    cfgs = dec2bin(0:15) == '1'; reps = unique([F.rep]);
    V = zeros(numel(reps), 16, numel(outs));
    for ir = 1:numel(reps)
        for ic = 1:16
            k = find([F.rep] == reps(ir) & [F.B] == cfgs(ic,1) & [F.A] == cfgs(ic,2) & [F.G] == cfgs(ic,3) & [F.S] == cfgs(ic,4));
            V(ir, ic, :) = [log(F(k).PRfm) log(F(k).IRRfm) 100*(F(k).effF - F(k).effM) F(k).Gabs];
        end
    end
    Cfg = struct([]);
    for ic = 1:16
        k = [F.B] == cfgs(ic,1) & [F.A] == cfgs(ic,2) & [F.G] == cfgs(ic,3) & [F.S] == cfgs(ic,4);
        p = ci([F(k).PRfm]); r = ci([F(k).IRRfm]); e = ci(100*([F(k).effF] - [F(k).effM]));
        Cfg = addrow(Cfg, 'B', cfgs(ic,1), 'A', cfgs(ic,2), 'G', cfgs(ic,3), 'S', cfgs(ic,4), ...
            'PR', p(1), 'PR_lo', p(2), 'PR_hi', p(3), 'IRR', r(1), 'IRR_lo', r(2), 'IRR_hi', r(3), ...
            'dEff', e(1), 'dEff_lo', e(2), 'dEff_hi', e(3), 'effF', mean([F(k).effF]), 'effM', mean([F(k).effM]), ...
            'effPreg', mean([F(k).effCovPreg]), 'Gabs', mean([F(k).Gabs]), 'PRadult', mean([F(k).PRadult]), ...
            'PresF', mean([F(k).PresF]), 'PresM', mean([F(k).PresM]), 'riskRec', mean([F(k).riskRecStd]), ...
            'riskNon', mean([F(k).riskNonStd]), 'RRrec', mean([F(k).RRrecStd]));
    end
    struct2csv('results/summary_A2_configurations.csv', Cfg);
    names = {'Behaviour (B)', 'Programme access (A)', 'Pregnancy biology (G)', 'Sex-specific clearance (S)', 'Interactions'};
    Sh = struct([]);
    for io = 1:numel(outs)
        phi = zeros(numel(reps), 4); inter = zeros(numel(reps), 1);
        for ir = 1:numel(reps)
            v = V(ir, :, io); phi(ir,:) = shapley(v, cfgs);
            main = arrayfun(@(j) v(bin2row(cfgs, (1:4) == j)) - v(1), 1:4);
            inter(ir) = (v(16) - v(1)) - sum(main);
        end
        tot = ci(V(:,16,io) - V(:,1,io));
        for j = 1:5
            if j <= 4, q = ci(phi(:,j)); else, q = ci(inter); end
            Sh = addrow(Sh, 'outcome', outs{io}, 'component', names{j}, 'mean', q(1), 'lo', q(2), 'hi', q(3), ...
                'total', tot(1), 'total_lo', tot(2), 'total_hi', tot(3));
            fprintf(fid, '   %-7s %-28s %+.4f (%+.4f,%+.4f)   [total %+.4f]\n', outs{io}, names{j}, q, tot(1));
        end
    end
    struct2csv('results/summary_A2_shapley.csv', Sh);
end

% ------------------------------------------------------------ A3 pregnancy sensitivity
if exist('results/A3_pregnancy.csv', 'file')
    G3 = csv2struct('results/A3_pregnancy.csv'); vs = unique({G3.variant}, 'stable');
    A3 = struct([]); fprintf(fid, '\nA3 pregnancy sensitivity (pregnant vs non-pregnant women 15-44 y)\n');
    for i = 1:numel(vs)
        k = strcmp({G3.variant}, vs{i});
        p = ci([G3(k).PRpreg]); r = ci([G3(k).IRRpreg]);
        A3 = addrow(A3, 'variant', vs{i}, 'PRpreg', p(1), 'PRpreg_lo', p(2), 'PRpreg_hi', p(3), ...
            'IRRpreg', r(1), 'IRRpreg_lo', r(2), 'IRRpreg_hi', r(3), 'rdtPreg', mean([G3(k).rdtPreg]), ...
            'rdtNonPreg', mean([G3(k).rdtNonPreg]), 'PRfm', mean([G3(k).PRfm]));
        fprintf(fid, '   %-14s PR %.2f (%.2f-%.2f)  IRR %.2f (%.2f-%.2f)  RDT preg %.3f non-preg %.3f\n', vs{i}, p, r, ...
            mean([G3(k).rdtPreg]), mean([G3(k).rdtNonPreg]));
    end
    struct2csv('results/summary_A3.csv', A3);
end

% ------------------------------------------------------------ A4 ITN use mode
if exist('results/A4_itn_use.csv', 'file')
    U4 = csv2struct('results/A4_itn_use.csv'); A4 = struct([]);
    fprintf(fid, '\nA4 ITN use (habitual vs nightly)\n');
    for nightly = [0 1]
        k = [U4.nightly] == nightly & strcmp({U4.rule}, 'policyProb');
        kP = [U4.nightly] == nightly & strcmp({U4.rule}, 'P'); kE = [U4.nightly] == nightly & strcmp({U4.rule}, 'E');
        o = ci([U4(k).ORu5]); rr = ci([U4(k).RRrecStd]); dEP = ci(100*([U4(kE).prev]./[U4(kP).prev] - 1));
        ir = ci([U4(k).IRRinf]);
        A4 = addrow(A4, 'nightly', nightly, 'ORu5', o(1), 'ORu5_lo', o(2), 'ORu5_hi', o(3), ...
            'riskRec', mean([U4(k).riskRecStd]), 'riskNon', mean([U4(k).riskNonStd]), ...
            'RRrec', rr(1), 'RRrec_lo', rr(2), 'RRrec_hi', rr(3), 'infRec', mean([U4(k).infRec]), ...
            'infNon', mean([U4(k).infNon]), 'IRRinf', ir(1), 'IRRinf_lo', ir(2), 'IRRinf_hi', ir(3), ...
            'EvsP', dEP(1), 'EvsP_lo', dEP(2), 'EvsP_hi', dEP(3));
        fprintf(fid, '   nightly=%d OR %.2f (%.2f-%.2f) RR recipients %.2f (%.2f-%.2f) risk rec %.3f non %.3f; infections/person rec %.2f non %.2f IRR %.2f (%.2f-%.2f); E vs P prevalence %+.1f%% (%+.1f,%+.1f)\n', ...
            nightly, o, rr, mean([U4(k).riskRecStd]), mean([U4(k).riskNonStd]), mean([U4(k).infRec]), mean([U4(k).infNon]), ir, dEP);
    end
    struct2csv('results/summary_A4.csv', A4);
end

% ------------------------------------------------------------ A5 allocation
if exist('results/A5_allocation.csv', 'file')
    L = csv2struct('results/A5_allocation.csv'); rules = {'P','U','D','E','A','H','T'};
    kP = strcmp({L.rule}, 'P'); A5 = struct([]);
    fprintf(fid, '\nA5 allocation (paired against P)\n');
    for i = 1:numel(rules)
        k = strcmp({L.rule}, rules{i});
        pv = ci(100*[L(k).prev]); in = ci([L(k).inc]);
        dp = ci(100*([L(k).prev]./[L(kP).prev] - 1)); di = ci(100*([L(k).inc]./[L(kP).inc] - 1));
        dz = ci(100*([L(k).Ztot]./[L(kP).Ztot] - 1)); du = ci(100*([L(k).incU5]./[L(kP).incU5] - 1));
        dpp = ci(100*([L(k).rdtPreg] - [L(kP).rdtPreg]));
        sev = 1000*sum([L(k).nSev])/sum([L(k).py]);
        A5 = addrow(A5, 'rule', rules{i}, 'prev', pv(1), 'inc', in(1), ...
            'dPrev', dp(1), 'dPrev_lo', dp(2), 'dPrev_hi', dp(3), 'dInc', di(1), 'dInc_lo', di(2), 'dInc_hi', di(3), ...
            'dZ', dz(1), 'dZ_lo', dz(2), 'dZ_hi', dz(3), 'dIncU5', du(1), 'dIncU5_lo', du(2), 'dIncU5_hi', du(3), ...
            'dPregRDT', dpp(1), 'dPregRDT_lo', dpp(2), 'dPregRDT_hi', dpp(3), 'sevInc', sev, ...
            'effF', mean([L(k).effF]), 'effM', mean([L(k).effM]), 'effPreg', mean([L(k).effCovPreg]), ...
            'Gabs', mean([L(k).Gabs]), 'TI', mean([L(k).TI]));
        fprintf(fid, '   %s prev %.2f%% inc %.1f | dPrev %+.1f%% (%+.1f,%+.1f) dInc %+.1f%% dZ %+.1f%% dIncU5 %+.1f%% dPregRDT %+.2f pp | sev %.2f eff F %.2f M %.2f preg %.2f Gabs %.1f\n', ...
            rules{i}, pv(1), in(1), dp, di(1), dz(1), du(1), dpp(1), sev, mean([L(k).effF]), mean([L(k).effM]), mean([L(k).effCovPreg]), mean([L(k).Gabs]));
    end
    struct2csv('results/summary_A5.csv', A5);
end

% ------------------------------------------------------------ A6 biting ratio
if exist('results/A6_biting_ratio.csv', 'file')
    W = csv2struct('results/A6_biting_ratio.csv'); A6 = struct([]);
    fprintf(fid, '\nA6 biting ratio\n');
    Rv = unique([W.Rout]);
    for iR = 1:numel(Rv)
        for md = {'fixed', 'recal'}
            if Rv(iR) == 1 && strcmp(md{1}, 'recal'), continue; end
            k = [W.Rout] == Rv(iR) & strcmp({W.mode}, md{1});
            if ~any(k), continue; end
            g = @(run) W(k & strcmp({W.run}, run));
            d0 = g('delta0.00'); d1 = g('delta-0.15'); d2 = g('delta-1.00'); rP = g('ruleP'); rE = g('ruleE');
            b1 = ci([d1.dPrev] - [d0.dPrev]); b2 = ci([d2.dPrev] - [d0.dPrev]);
            ep = ci(100*([rE.prev]./[rP.prev] - 1));
            A6 = addrow(A6, 'Rout', Rv(iR), 'mode', md{1}, 'm', d0(1).m, 'prev', mean([rP.prev]), 'u5rdt', mean([rP.u5rdt]), ...
                'outF', mean([d2.outF]), 'outM', mean([d2.outM]), 'beh015', b1(1), 'beh015_lo', b1(2), 'beh015_hi', b1(3), ...
                'beh100', b2(1), 'beh100_lo', b2(2), 'beh100_hi', b2(3), 'EvsP', ep(1), 'EvsP_lo', ep(2), 'EvsP_hi', ep(3));
            fprintf(fid, '   R %.2f %-5s m %.1f prev %.3f u5 %.3f | beh(-0.15) %+.2f (%+.2f,%+.2f) beh(-1) %+.2f (%+.2f,%+.2f) | E vs P %+.1f%% (%+.1f,%+.1f)\n', ...
                Rv(iR), md{1}, d0(1).m, mean([rP.prev]), mean([rP.u5rdt]), b1, b2, ep);
        end
    end
    struct2csv('results/summary_A6.csv', A6);
end

% ------------------------------------------------------------ A7 uncertainty
merge_parts('results/A7_uncertainty');          % joins A7 runs split across sessions
if exist('results/A7_uncertainty.csv', 'file')
    H = csv2struct('results/A7_uncertainty.csv'); n = numel(H);
    pr = @(x) mean(x);
    A7 = struct('n', n, ...
        'p_E_beats_P', pr([H.prev_E] < [H.prev_P]), 'p_A_beats_P', pr([H.prev_A] < [H.prev_P]), ...
        'p_U_beats_P', pr([H.prev_U] < [H.prev_P]), 'p_D_beats_P', pr([H.prev_D] < [H.prev_P]), ...
        'p_E_best', pr(([H.prev_E] <= min([[H.prev_P];[H.prev_U];[H.prev_D];[H.prev_A]]))), ...
        'p_clear_lowers_F', pr([H.PR_clear] < 1), 'p_occ_lowers_F', pr([H.PR_occ] < 1), ...
        'p_clear_gt_occ', pr(abs(log([H.PR_clear])) > abs(log([H.PR_occ]))), ...
        'p_Ppreg_better_than_E', pr([H.preg_P] < [H.preg_E]), 'p_E_raises_U5', pr([H.incU5_E] > [H.incU5_P]), ...
        'med_EvsP', median(100*([H.prev_E]./[H.prev_P] - 1)), 'med_AvsP', median(100*([H.prev_A]./[H.prev_P] - 1)), ...
        'med_PRclear', median([H.PR_clear]), 'med_PRocc', median([H.PR_occ]));
    fprintf(fid, '\nA7 uncertainty (%d draws)\n', n);
    f = fieldnames(A7); for i = 1:numel(f), fprintf(fid, '   %-22s %.3f\n', f{i}, A7.(f{i})); end
    struct2csv('results/summary_A7.csv', A7);
end

% ------------------------------------------------------------ A9 mechanism checks
if exist('results/A9_immunity_runs.csv', 'file')
    H = csv2struct('results/A9_immunity_runs.csv'); A9 = struct([]);
    imm = [H.immType]; cfg = {H.config};
    for it = unique(imm)
        for c = {'equality', 'clearance only', 'reference'}
            k = imm == it & strcmp(cfg, c{1}); row = struct('immType', it, 'config', c{1});
            for f = {'PRfm', 'IRRfm', 'prev', 'inc', 'suscF', 'suscM', 'durIAF', 'durIAM', ...
                     'infPYF', 'infPYM', 'clinPYF', 'clinPYM', 'pCF', 'pCM'}
                q = ci([H(k).(f{1})]); row.(f{1}) = q(1); row.([f{1} '_lo']) = q(2); row.([f{1} '_hi']) = q(3);
            end
            q = ci([H(k).infPYF]./[H(k).infPYM]);   row.infRatio = q(1);  row.infRatio_lo = q(2);  row.infRatio_hi = q(3);
            q = ci([H(k).clinPYF]./[H(k).clinPYM]); row.clinRatio = q(1); row.clinRatio_lo = q(2); row.clinRatio_hi = q(3);
            q = ci([H(k).durIAF]./[H(k).durIAM]);   row.durRatio = q(1);  row.durRatio_lo = q(2);  row.durRatio_hi = q(3);
            q = ci([H(k).pCF]./[H(k).pCM]);         row.pCRatio = q(1);   row.pCRatio_lo = q(2);   row.pCRatio_hi = q(3);
            if isempty(A9), A9 = row; else, A9(end+1) = orderfields(row, A9(1)); end %#ok<AGROW>
        end
    end
    struct2csv('results/summary_A9_immunity.csv', A9);
end
if exist('results/A9_itn_components.csv', 'file')
    H = csv2struct('results/A9_itn_components.csv'); Q = struct();
    for f = {'avertedUser', 'infHold', 'infNonHold', 'infNoITN_hold', 'infNoITN_non', 'RRdirect', ...
             'RRindirect', 'RRtotal', 'prevITN', 'prevNoITN', 'incITN', 'incNoITN', 'u5share', 'pregShare'}
        q = ci([H.(f{1})]); Q.(f{1}) = q(1); Q.([f{1} '_lo']) = q(2); Q.([f{1} '_hi']) = q(3);
    end
    q = ci(1 - [H.prevITN]./[H.prevNoITN]); Q.prevReduction = q(1); Q.prevReduction_lo = q(2); Q.prevReduction_hi = q(3);
    q = ci(1 - [H.incITN]./[H.incNoITN]);   Q.incReduction = q(1);  Q.incReduction_lo = q(2);  Q.incReduction_hi = q(3);
    struct2csv('results/summary_A9_itn.csv', Q);
end

end

% =========================================================================
function q = ci(v)
v = v(:); v = v(isfinite(v)); n = numel(v);
if n < 2, q = [mean(v) NaN NaN]; return; end
h = tq(n-1)*std(v)/sqrt(n); q = [mean(v) mean(v)-h mean(v)+h];
end

function t = tq(df)
% 97.5th percentile of Student's t (table for small df, normal limit otherwise)
tab = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 2.201 2.179 2.160 2.145 2.131 ...
       2.120 2.110 2.101 2.093 2.086 2.080 2.074 2.069 2.064 2.060 2.056 2.052 2.048 2.045 2.042];
if df <= numel(tab), t = tab(df); else, t = 1.96 + 2.4/df; end
end

function i = bin2row(cfgs, b)
i = find(all(cfgs == repmat(b, size(cfgs,1), 1), 2));
end

function phi = shapley(v, cfgs)
% Shapley value of each of 4 factors for outcome v(config) relative to v(0000)
k = 4; phi = zeros(1,k);
for j = 1:k
    for ic = 1:size(cfgs,1)
        s = cfgs(ic,:);
        if s(j), continue; end
        sz = sum(s); wgt = factorial(sz)*factorial(k - sz - 1)/factorial(k);
        t = s; t(j) = 1;
        phi(j) = phi(j) + wgt*(v(bin2row(cfgs, t)) - v(ic));
    end
end
end

function merge_parts(stem)
% join files stem_<a>_<b>.csv written by parallel partial runs into stem.csv
d = dir([stem '_*_*.csv']);
if isempty(d), return; end
[p, ~] = fileparts(stem); T = struct([]);
names = sort({d.name});
for i = 1:numel(names)
    H = csv2struct(fullfile(p, names{i}));
    if isempty(T), T = H; else, T = [T(:); H(:)]'; end %#ok<AGROW>
end
struct2csv([stem '.csv'], T);
end
