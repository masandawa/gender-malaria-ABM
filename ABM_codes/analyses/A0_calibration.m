function m_star = A0_calibration(X)
%A0_CALIBRATION  Single-target calibration of the mosquito-to-human ratio m
%   and a profile analysis of parameter identifiability.
%
%   Target: RDT-detectable prevalence of 6% in children aged 6-59 months,
%   Morogoro region, TDHS-MIS 2022 (acceptance 4-8%). Morogoro includes the
%   Kilombero valley where the entomological and behavioural inputs were
%   measured. Objective: squared error on the log scale. m* is obtained by
%   inverting a log-log regression of simulated prevalence on m; its 95%
%   interval comes from 1000 bootstrap resamples of the realisations.
%   The profile repeats the calibration with alternative values of parameters
%   that a single prevalence target cannot identify.
P0 = base_params(X);
runs = struct([]);
for im = 1:numel(X.calGrid)
    for k = 1:X.calReps
        R = cal_run(X, X.calGrid(im), struct(), 1000*im + k);
        runs = addrow(runs, 'set', 'reference', 'm', X.calGrid(im), 'rep', k, R);
    end
    fprintf('[A0] m = %g done\n', X.calGrid(im));
end
[m_star, b] = fit_m([runs.m], [runs.u5rdt], P0.calTarget);
mLo = inv_fit(b, P0.calLo); mHi = inv_fit(b, P0.calHi);
% bootstrap interval for m*
mm = [runs.m]; yy = [runs.u5rdt]; B = 1000; mb = NaN(B,1);
set_seed(77);
for i = 1:B
    id = randi(numel(mm), numel(mm), 1);
    mb(i) = fit_m(mm(id), yy(id), P0.calTarget);
end
mb = sort(mb(isfinite(mb)));
ci = mb(max(1, round([0.025 0.975]*numel(mb))));
inAcc = arrayfun(@(m) mean(yy(mm == m) >= P0.calLo & yy(mm == m) <= P0.calHi), unique(mm));
S = struct('m_star', m_star, 'boot_lo', ci(1), 'boot_hi', ci(2), 'acc_lo', mLo, 'acc_hi', mHi, ...
           'slope', b(2), 'target', P0.calTarget);
fprintf('[A0] m* = %.2f (bootstrap 95%% %.2f-%.2f; acceptance %.2f-%.2f)\n', m_star, ci(1), ci(2), mLo, mHi);

% ---------------------------------------------------------- profile analysis
alt = {'d_A = 150 d', 'dA', 150; 'd_A = 250 d', 'dA', 250; ...
       'eta_N = 0.50', 'etaN', 0.50; 'eta_N = 0.90', 'etaN', 0.90; ...
       'b_vh x 0.75', 'bvh', 0.75*P0.bvh};
prof = struct([]);
for ia = 1:size(alt,1)
    ov = struct(alt{ia,2}, alt{ia,3});
    grid = m_star*[0.7 1.0 1.4];
    for im = 1:3
        for k = 1:X.profReps
            R = cal_run(X, grid(im), ov, 50000 + 100*ia + 10*im + k);
            runs = addrow(runs, 'set', alt{ia,1}, 'm', grid(im), 'rep', k, R);
        end
    end
    kk = strcmp({runs.set}, alt{ia,1});
    [ma, ba] = fit_m([runs(kk).m], [runs(kk).u5rdt], P0.calTarget);
    % implied outputs at the profile optimum (log-linear interpolation)
    inc = interp_out([runs(kk).m], [runs(kk).inc], ma);
    prv = interp_out([runs(kk).m], [runs(kk).prev], ma);
    prof = addrow(prof, 'set', alt{ia,1}, 'm_star', ma, 'slope', ba(2), 'inc', inc, 'prev', prv);
    fprintf('[A0] profile %s: m* = %.2f\n', alt{ia,1}, ma);
end
kk = strcmp({runs.set}, 'reference');
prof = addrow(prof, 'set', 'reference', 'm_star', m_star, 'slope', b(2), ...
    'inc', interp_out([runs(kk).m], [runs(kk).inc], m_star), ...
    'prev', interp_out([runs(kk).m], [runs(kk).prev], m_star));
struct2csv('results/A0_calibration_runs.csv', runs);
struct2csv('results/A0_calibration_summary.csv', S);
struct2csv('results/A0_profile.csv', prof);
acc = struct('m', num2cell(unique(mm)), 'p_accept', num2cell(inAcc));
struct2csv('results/A0_acceptance.csv', acc);
end

% =========================================================================
function R = cal_run(X, m, ov, seed)
P = base_params(X, m); [P, C] = scenario(P, 'REF');
f = fieldnames(ov); for j = 1:numel(f), P.(f{j}) = ov.(f{j}); end
set_seed(seed); S = abm_init(P);
S = abm_run(P, S, X.burn, C);
[~, O] = abm_run(P, S, 365, C);
R = summarize_run(P, O);
R = struct('u5rdt', R.u5rdt, 'rdt', R.rdt, 'prev', R.prev, 'inc', R.inc);
end

function [m, b] = fit_m(mv, y, target)
ok = y > 0;
if numel(unique(mv(ok))) < 2, m = NaN; b = [NaN; NaN]; return; end
Xd = [ones(sum(ok),1) log(mv(ok)')]; b = Xd \ log(y(ok)');
m = inv_fit(b, target);
end

function m = inv_fit(b, y)
m = exp((log(y) - b(1))/b(2));
end

function v = interp_out(mv, y, m)
ok = y > 0;
b = [ones(sum(ok),1) log(mv(ok)')] \ log(y(ok)');
v = exp(b(1) + b(2)*log(m));
end
