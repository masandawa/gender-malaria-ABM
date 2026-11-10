function make_figures()
%MAKE_FIGURES  Figures 2-6 of the manuscript from results/*.csv.
%   Figures are drawn at their printed size (17 cm wide) with 7.5-9 pt text and
%   saved as 600-dpi PNG and vector PDF in figures/. Category labels are placed
%   on the vertical axis so that no label needs rotating.
if exist('OCTAVE_VERSION', 'builtin'), graphics_toolkit('gnuplot'); end
if ~exist('figures', 'dir'), mkdir('figures'); end
fig2_calibration();
fig3_behaviour();
fig4_mechanisms();
fig5_itn();
fig6_robustness();
end

% =========================================================================
function fig2_calibration()
R  = csv2struct('results/A0_calibration_runs.csv');
S  = csv2struct('results/A0_calibration_summary.csv');
Pf = csv2struct('results/A0_profile.csv');
E  = csv2struct('results/summary_evaluation.csv');
P  = abm_params();
f = fig_new(17, 15);
% (a) calibration curve
ax = panel([0.09 0.62 0.36 0.31], 'a', 'Single-target calibration');
k = strcmp({R.set}, 'reference'); m = [R(k).m]; y = 100*[R(k).u5rdt];
xr = [min(m) - 1, max(m) + 1]; ytop = max([y 100*P.calHi])*1.35;
fill([xr(1) xr(2) xr(2) xr(1)], 100*[P.calLo P.calLo P.calHi P.calHi], okabe('band'), 'EdgeColor', 'none');
plot(m, y, 'o', 'Color', okabe('male'), 'MarkerSize', 3.5, 'MarkerFaceColor', okabe('sky'));
mm = linspace(xr(1), xr(2), 60);
plot(mm, 100*P.calTarget*(mm/S.m_star).^S.slope, 'k-', 'LineWidth', 1.3);
plot([S.m_star S.m_star], [0 ytop], '--', 'Color', okabe('female'), 'LineWidth', 1.1);
text(min(S.m_star, xr(2) - 3) + 0.3, 0.93*ytop, sprintf('m^* = %.1f', S.m_star), 'FontSize', 8.5);
xlim(xr); ylim([0 ytop]);
xlabel('Mosquito-to-human ratio m'); ylabel('RDT prevalence, 6-59 months (%)');
mlegend(ax, {'Target 4-8%', 'Realisation', 'Fitted'}, ...
    {{okabe('band'), 'bar', 'none'}, {okabe('sky'), 'o', 'none'}, {[0 0 0], 'none', '-'}}, [0.03 0.95], 0.09);
% (b) profile
ax = panel([0.70 0.62 0.27 0.31], 'b', 'Profile: alternative fits');
o = [find(strcmp({Pf.set}, 'reference')) find(~strcmp({Pf.set}, 'reference'))]; Pf = Pf(o); n = numel(Pf);
lab = {Pf.set}; lab = strrep(lab, 'reference', 'Reference'); lab = strrep(lab, 'eta_N', '\eta_N');
lab = strrep(lab, 'b_vh x 0.75', 'b_{vh} \times 0.75'); lab = strrep(lab, ' = ', ' ');
b = barh(n:-1:1, [Pf.m_star], 0.6); set(b, 'FaceColor', [0.78 0.85 0.92], 'EdgeColor', [0.3 0.3 0.3]);
for i = 1:n
    text(Pf(i).m_star + 0.6, n - i + 1, sprintf('%.0f', Pf(i).inc), 'FontSize', 7.5, 'VerticalAlignment', 'middle');
end
set(ax, 'YTick', 1:n, 'YTickLabel', fliplr(lab)); ylim([0.4 n + 0.6]); xlim([0 1.3*max([Pf.m_star])]);
xlabel('Calibrated m^* (label: incidence)');
% (c) proportions
ax = panel([0.30 0.10 0.20 0.38], 'c', 'Proportions');
Vp = {'u5rdt', 0.06, 0.04, 0.08, 'RDT prev. 6-59 m (C)'; 'piI', 0.652, 0.652, 0.754, 'Indoor exposure share'; ...
      'use', 0.59, 0.59, 0.59, 'ITN use, all (I)'; 'effCovPreg', 0.65, 0.64, 0.66, 'ITN use, pregnant (I)'};
forest(ax, E, Vp); xlim([0 0.9]); xlabel('Proportion');
mlegend(ax, {'Observed', 'Model'}, {{[0 0.35 0.25], 'd', '-'}, {[0 0 0], 'o', '-'}}, [0.40 0.98], 0.075);
% (d) ratios
ax = panel([0.77 0.10 0.20 0.38], 'd', 'Ratios (log scale)');
Vr = {'PRfmRDT', 0.805, 0.763, 0.849, 'F:M prevalence, all ages'; 'PRu5', 1.00, 1.00, 1.00, 'F:M prevalence, <5 y'; ...
      'PregU5', 1.16, 0.90, 1.51, 'Pregnant : children (R)'; 'ORu5', 2.73, 1.51, 4.91, 'OR, no ITN, <5 y (R)'};
forest(ax, E, Vr); set(ax, 'XScale', 'log'); xlim([0.45 6]);
set(ax, 'XTick', [0.5 1 2 4], 'XTickLabel', {'0.5', '1', '2', '4'}); xlabel('Ratio');
fig_save(f, 'Figure2');
end

function forest(ax, E, V)
n = size(V,1);
for i = 1:n
    yv = n - i + 1; q = E(strcmp({E.quantity}, V{i,1}));
    plot(ax, [V{i,3} V{i,4}], [yv yv] + 0.17, '-', 'Color', okabe('green'), 'LineWidth', 3.5);
    plot(ax, V{i,2}, yv + 0.17, 'd', 'Color', [0 0.35 0.25], 'MarkerFaceColor', [0 0.35 0.25], 'MarkerSize', 4.5);
    plot(ax, [q.lo q.hi], [yv yv] - 0.17, 'k-', 'LineWidth', 1.4);
    plot(ax, q.mean, yv - 0.17, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3.5);
end
set(ax, 'YTick', 1:n, 'YTickLabel', fliplr(V(:,5)')); ylim(ax, [0.4 n + 0.75]);
end

% =========================================================================
function fig3_behaviour()
P  = abm_params();
A1 = csv2struct('results/summary_A1.csv');
Cs = csv2struct('results/summary_A1_concentration.csv');
f = fig_new(17, 13.5);
h = 1:14; lab = {'18','','20','','22','','00','','02','','04','','06',''};
sg = @(x) 1./(1 + exp(-x)); lA = log(P.poutA./(1 - P.poutA));
% (a) activity
ax = panel([0.08 0.61 0.38 0.31], 'a', 'Probability of being outdoors');
plot(h, sg(lA + P.deltaA/2), '-o', 'Color', okabe('female'), 'MarkerSize', 3, 'MarkerFaceColor', okabe('female'));
plot(h, sg(lA - P.deltaA/2), '-s', 'Color', okabe('male'), 'MarkerSize', 3, 'MarkerFaceColor', okabe('male'));
plot(h, sg(lA + 0.5), '--', 'Color', okabe('male'));
plot(h, P.poutY, '-^', 'Color', okabe('grey'), 'MarkerSize', 3);
set(ax, 'XTick', h, 'XTickLabel', lab); xlim([0.5 14.5]); ylim([0 1.05]);
xlabel('Hour starting'); ylabel('Probability');
mlegend(ax, {'Women, \delta_A = -0.15', 'Men, \delta_A = -0.15', 'Men, \delta_A = -1.0', 'Children <5 y'}, ...
    {{okabe('female'), 'o', '-'}, {okabe('male'), 's', '-'}, {okabe('male'), 'none', '--'}, {okabe('grey'), '^', '-'}}, ...
    [0.33 0.80], 0.085);
% (b) biting
ax = panel([0.58 0.61 0.38 0.31], 'b', 'Relative mosquito biting');
b = bar(h, 100*P.w, 0.7); set(b, 'FaceColor', [0.70 0.70 0.70], 'EdgeColor', [0.3 0.3 0.3]);
set(ax, 'XTick', h, 'XTickLabel', lab); xlim([0.5 14.5]); ylim([0 1.3*max(100*P.w)]);
xlabel('Hour starting'); ylabel('Share of nightly bites (%)');
text(1, 1.2*max(100*P.w), sprintf('Reference R_{out/in} = %.2f', P.Rout), 'FontSize', 8);
% (c) behavioural component
ax = panel([0.10 0.10 0.36 0.31], 'c', 'Behavioural female-male difference');
d = [A1.delta];
plot([-1.1 0.4], [0 0], ':', 'Color', okabe('grey'));
eb(d, [A1.dBeh], [A1.dBeh_lo], [A1.dBeh_hi], 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3.5, 'LineWidth', 1.1);
xlim([-1.1 0.4]); xlabel('\delta_A (negative: men more often outdoors)'); ylabel('Change in F-M prevalence (pp)');
% (d) concentration
ax = panel([0.58 0.10 0.38 0.31], 'd', 'Concentration of events');
Y = 100*[Cs.qB50 Cs.qT50; Cs.qB80 Cs.qT80];
b = bar(1:2, Y, 0.8); set(b(1), 'FaceColor', okabe('orange')); set(b(2), 'FaceColor', okabe('sky'));
set(ax, 'XTick', 1:2, 'XTickLabel', {'50% of events', '80% of events'}); xlim([0.4 2.6]);
ylabel('Smallest population fraction (%)'); ylim([0 1.45*max(Y(:))]);
mlegend(ax, {'Clinical episodes', 'Attributed transmission', sprintf('Overlap J_+ = %.2f', Cs.J)}, ...
    {{okabe('orange'), 'bar', 'none'}, {okabe('sky'), 'bar', 'none'}, {[1 1 1], 'none', 'none'}}, [0.03 0.97], 0.085);
fig_save(f, 'Figure3');
end

% =========================================================================
function fig4_mechanisms()
Sh = csv2struct('results/summary_A2_shapley.csv');
A3 = csv2struct('results/summary_A3.csv');
f = fig_new(17, 13.5);
comps = {'Behaviour (B)', 'Programme access (A)', 'Pregnancy biology (G)', 'Sex-specific clearance (S)', 'Interactions'};
short = {'Behaviour', 'Programme', 'Pregnancy', 'Clearance', 'Interaction'};
outs = {'logPR', 'Prevalence ratio, log(F:M)'; 'logIRR', 'Incidence ratio, log(F:M)'; 'dEff', 'Effective ITN coverage, F-M'};
xl = {'Shapley contribution', 'Shapley contribution', 'Shapley contribution (pp)'};
pos = {[0.15 0.61 0.31 0.31], [0.65 0.61 0.31 0.31], [0.15 0.10 0.31 0.31]};
let = 'abc';
for io = 1:3
    ax = panel(pos{io}, let(io), outs{io,2});
    s = Sh(strcmp({Sh.outcome}, outs{io,1}));
    for j = 1:5
        q = s(strcmp({s.component}, comps{j}));
        col = okabe('male'); if q.mean > 0, col = okabe('female'); end
        b = barh(6 - j, q.mean, 0.6); set(b, 'FaceColor', col, 'EdgeColor', 'none');
        plot([q.lo q.hi], [6 - j 6 - j], 'k-', 'LineWidth', 1.1);
    end
    yl = [0.4 5.6]; plot([0 0], yl, 'k-', 'LineWidth', 0.6);
    set(ax, 'YTick', 1:5, 'YTickLabel', fliplr(short)); ylim(yl);
    xlabel(sprintf('%s (total %+.2f)', xl{io}, s(1).total));
end
% (d) pregnancy sensitivity
ax = panel([0.65 0.10 0.31 0.31], 'd', 'Pregnant vs non-pregnant women');
n = numel(A3); yv = n:-1:1;
for i = 1:n
    plot([A3(i).PRpreg_lo A3(i).PRpreg_hi], [yv(i) yv(i)] + 0.17, '-', 'Color', okabe('preg'), 'LineWidth', 1.2);
    plot([A3(i).IRRpreg_lo A3(i).IRRpreg_hi], [yv(i) yv(i)] - 0.17, '-', 'Color', 'k', 'LineWidth', 1.2);
end
plot([A3.PRpreg], yv + 0.17, 'o', 'Color', okabe('preg'), 'MarkerFaceColor', okabe('preg'), 'MarkerSize', 3.5);
plot([A3.IRRpreg], yv - 0.17, 's', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3.5);
hi = max([[A3.PRpreg_hi] [A3.IRRpreg_hi]]); xlim([0 1.75*hi]);
plot([1 1], [0.4 n + 1.6], ':', 'Color', okabe('grey'));
lab = {A3.variant}; lab = strrep(lab, 'reference', 'Reference'); lab = strrep(lab, 'kappa_P', '\kappa_P');
lab = strrep(lab, 'r_C', 'r^P_C'); lab = strrep(lab, 'IPTp', 'IPTp eff.');
set(ax, 'YTick', 1:n, 'YTickLabel', fliplr(lab)); ylim([0.4 n + 1.6]);
xlabel('Ratio, pregnant : non-pregnant');
mlegend(ax, {'RDT prevalence', 'Clinical incidence'}, {{okabe('preg'), 'o', '-'}, {[0 0 0], 's', '-'}}, [0.62 0.98], 0.07);
fig_save(f, 'Figure4');
end

% =========================================================================
function fig5_itn()
C  = csv2struct('results/summary_A2_configurations.csv');
A4 = csv2struct('results/summary_A4.csv');
A5 = csv2struct('results/summary_A5.csv');
f = fig_new(17, 14);
% (a) effective coverage
ax = panel([0.08 0.62 0.38 0.30], 'a', 'Realised effective ITN coverage');
sel = [0 0 0 0; 1 0 0 0; 0 1 0 0; 1 1 0 0; 1 1 1 1]; names = {'Equality', '+B', '+A', '+B+A', 'Reference'};
Y = zeros(5,3);
for i = 1:5
    k = [C.B] == sel(i,1) & [C.A] == sel(i,2) & [C.G] == sel(i,3) & [C.S] == sel(i,4);
    Y(i,:) = 100*[C(k).effF C(k).effM C(k).effPreg];
end
b = bar(1:5, Y, 0.85);
set(b(1), 'FaceColor', okabe('female')); set(b(2), 'FaceColor', okabe('male')); set(b(3), 'FaceColor', okabe('preg'));
set(ax, 'XTick', 1:5, 'XTickLabel', names); xlim([0.4 5.6]); ylim([0 100]); ylabel('Person-nights under an ITN (%)');
mlegend(ax, {'Women', 'Men', 'Pregnant women'}, ...
    {{okabe('female'), 'bar', 'none'}, {okabe('male'), 'bar', 'none'}, {okabe('preg'), 'bar', 'none'}}, [0.03 0.97], 0.085);
% (b) standardised infection rate
ax = panel([0.58 0.62 0.38 0.30], 'b', 'Infections per person in 5 years');
Y = [A4(1).infRec A4(1).infNon; A4(2).infRec A4(2).infNon];
b = bar(1:2, Y, 0.7); set(b(1), 'FaceColor', okabe('green')); set(b(2), 'FaceColor', okabe('grey'));
set(ax, 'XTick', 1:2, 'XTickLabel', {'Habitual use (ref.)', 'Nightly use'}); xlim([0.4 2.6]); ylim([0 1.45*max(Y(:))]);
ylabel('Mean infections (standardised)');
for i = 1:2
    text(i, max(Y(i,:)) + 0.06*max(Y(:)), sprintf('RR %.2f', A4(i).IRRinf), 'HorizontalAlignment', 'center', 'FontSize', 8);
end
mlegend(ax, {'ITN holders', 'Non-holders'}, {{okabe('green'), 'bar', 'none'}, {okabe('grey'), 'bar', 'none'}}, [0.03 0.97], 0.085);
% (c) population outcomes vs policy-aligned
rules = {'U', 'D', 'E', 'A', 'H', 'T'};
rn = {'Uniform', 'Demographic', 'Exposure (ideal)', 'Activity proxy', 'Infection history', 'Transmission (ideal)'};
ax = panel([0.19 0.10 0.28 0.36], 'c', 'Change vs policy-aligned');
n = numel(rules); fl = {'dPrev', 'dInc'}; cols = {okabe('male'), okabe('female')}; mk = {'o', 's'};
for j = 1:2
    for i = 1:n
        a = A5(strcmp({A5.rule}, rules{i})); yv = n - i + 1 + (1.5 - j)*0.3;
        plot([a.([fl{j} '_lo']) a.([fl{j} '_hi'])], [yv yv], '-', 'Color', cols{j}, 'LineWidth', 1.2);
        plot(a.(fl{j}), yv, mk{j}, 'Color', cols{j}, 'MarkerFaceColor', cols{j}, 'MarkerSize', 3.5);
    end
end
vals = [[A5.dPrev_lo] [A5.dPrev_hi] [A5.dInc_lo] [A5.dInc_hi]]; lim = max(abs(vals))*1.15;
plot([0 0], [0.4 n + 1.4], 'k-', 'LineWidth', 0.6);
set(ax, 'YTick', 1:n, 'YTickLabel', fliplr(rn)); ylim([0.4 n + 1.4]); xlim([-lim lim]);
xlabel('Relative change (%)');
mlegend(ax, {'Prevalence', 'Clinical incidence'}, {{cols{1}, 'o', '-'}, {cols{2}, 's', '-'}}, [0.58 0.98], 0.07);
% (d) protection of pregnant women
ax = panel([0.69 0.10 0.27 0.36], 'd', 'Pregnant women: ITN coverage');
rr = [{'P'} rules]; rnn = [{'Policy-aligned'} rn]; m = numel(rr);
y = arrayfun(@(i) 100*A5(strcmp({A5.rule}, rr{i})).effPreg, 1:m);
b = barh(m:-1:1, y, 0.6); set(b, 'FaceColor', okabe('preg'), 'EdgeColor', [0.3 0.3 0.3]);
for i = 1:m, text(y(i) + 2, m - i + 1, sprintf('%.0f', y(i)), 'FontSize', 7.5, 'VerticalAlignment', 'middle'); end
set(ax, 'YTick', 1:m, 'YTickLabel', fliplr(rnn)); ylim([0.4 m + 0.6]); xlim([0 110]);
xlabel('Person-nights under an ITN (%)');
fig_save(f, 'Figure5');
end

% =========================================================================
function fig6_robustness()
A6 = csv2struct('results/summary_A6.csv');
A7 = csv2struct('results/summary_A7.csv');
H  = csv2struct('results/A7_uncertainty.csv');
f = fig_new(17, 14);
modes = {'fixed', 'recal'}; mn = {'Fixed density', 'Re-calibrated'}; ls = {'-', '--'}; mk = {'o', 's'};
cols = {okabe('male'), okabe('female')};
xt = [0.7 1 1.43 2.5 4.37]; xtl = {'0.70', '1.00', '1.43', '2.50', '4.37'};
fld = {'beh100', 'EvsP'}; ttl = {'Occupational contrast (\delta_A = -1.0)', 'Exposure (ideal) vs policy-aligned'};
yl = {'Change in F-M prevalence (pp)', 'Relative change in prevalence (%)'};
pos = {[0.08 0.62 0.38 0.30], [0.58 0.62 0.38 0.30]}; let = 'ab';
for p = 1:2
    ax = panel(pos{p}, let(p), ttl{p});
    for j = 1:2
        k = strcmp({A6.mode}, modes{j}) | [A6.Rout] == 1; a = A6(k); [~, o] = sort([a.Rout]); a = a(o);
        x = [a.Rout]*(1 + (j - 1.5)*0.05);
        for i = 1:numel(a), plot([x(i) x(i)], [a(i).([fld{p} '_lo']) a(i).([fld{p} '_hi'])], '-', 'Color', cols{j}); end
        plot(x, [a.(fld{p})], [ls{j} mk{j}], 'Color', cols{j}, 'MarkerFaceColor', cols{j}, 'MarkerSize', 3.5);
    end
    plot([0.6 5], [0 0], 'k:');
    set(ax, 'XScale', 'log', 'XTick', xt, 'XTickLabel', xtl); xlim([0.6 5]);
    xlabel('Outdoor:indoor biting ratio R_{out/in}'); ylabel(yl{p});
    v = [[A6.([fld{p} '_lo'])] [A6.([fld{p} '_hi'])] 0]; rg = max(v) - min(v);
    ylim([min(v) - 0.35*rg, max(v) + 0.1*rg]);
    mlegend(ax, mn, {{cols{1}, mk{1}, ls{1}}, {cols{2}, mk{2}, ls{2}}}, [0.03 0.17], 0.085);
end
% (c) probability that statements hold
ax = panel([0.30 0.10 0.17 0.36], 'c', sprintf('Joint uncertainty (%d draws)', A7.n));
pv = [A7.p_E_beats_P A7.p_E_raises_U5 A7.p_A_beats_P A7.p_D_beats_P A7.p_clear_lowers_F A7.p_occ_lowers_F A7.p_clear_gt_occ];
lb = {'Exposure (ideal) < policy', 'Exposure raises <5 y incidence', 'Activity proxy < policy', 'Demographic < policy', ...
      'Clearance only: F < M', 'Occupation only: F < M', '|Clearance| > |occupation|'};
b = barh(numel(pv):-1:1, 100*pv, 0.6); set(b, 'FaceColor', [0.70 0.78 0.88], 'EdgeColor', [0.3 0.3 0.3]);
for i = 1:numel(pv), text(100*pv(i) + 3, numel(pv) - i + 1, sprintf('%.0f%%', 100*pv(i)), 'FontSize', 7.5, 'VerticalAlignment', 'middle'); end
set(ax, 'YTick', 1:numel(pv), 'YTickLabel', fliplr(lb), 'XTick', 0:25:100); xlim([0 125]); ylim([0.4 numel(pv) + 0.6]);
xlabel('Draws (%)');
% (d) scatter
ax = panel([0.63 0.10 0.33 0.36], 'd', 'Female:male prevalence ratio by draw');
lims = [min([[H.PR_clear] [H.PR_occ]]) - 0.03, max([[H.PR_clear] [H.PR_occ] 1]) + 0.03];
plot(lims, lims, 'k:'); plot(lims, [1 1], ':', 'Color', okabe('grey')); plot([1 1], lims, ':', 'Color', okabe('grey'));
plot([H.PR_occ], [H.PR_clear], 'o', 'Color', okabe('male'), 'MarkerFaceColor', okabe('sky'), 'MarkerSize', 3.5);
xlim(lims); ylim(lims); xlabel('Occupational activity only'); ylabel('Sex-specific clearance only');
fig_save(f, 'Figure6');
end
