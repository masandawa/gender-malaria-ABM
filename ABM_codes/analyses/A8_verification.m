function T = A8_verification(X, m)
%A8_VERIFICATION  Internal verification of the implementation.
%   Each test runs a short simulation and checks an expected property.
%   Verification checks the code against its specification; it is not
%   evidence of epidemiological validity.
T = {};
P0 = base_params(X, m); [P0, C0] = scenario(P0, 'REF');
set_seed(1); S0 = abm_init(P0); S0 = abm_run(P0, S0, 365, C0);   % endemic start

% 1. b_vh = 0: no new human infections
P = P0; P.bvh = 0; [~, O] = abm_run(P, S0, 365, C0);
T(end+1,:) = {'b_vh = 0: no new human infections', sum(O.G(:,6)) == 0, sprintf('%d new infections', sum(O.G(:,6)))};
% 2. b_hv = 0: no new mosquito infections
P = P0; P.bhv = 0; [~, O] = abm_run(P, S0, 365, C0);
T(end+1,:) = {'b_hv = 0: no new mosquito infections', O.Ztot == 0, sprintf('%d mosquito infections', O.Ztot)};
% 3. population conservation and valid states
[S, O] = abm_run(P0, S0, 365, C0);
ok = numel(S.x) == P0.N && all(ismember(S.x, 0:4)) && S.Sv >= 0 && S.Ev >= 0 && S.Iv >= 0;
T(end+1,:) = {'Constant population size; valid human and mosquito states', ok, sprintf('N = %d; Sv,Ev,Iv = %d,%d,%d', numel(S.x), S.Sv, S.Ev, S.Iv)};
% 4. no infection seeding: no humans or mosquitoes infected -> no infection
S = S0; S.x(:) = 0; S.lat(:) = 0; S.Sv = S.Sv + S.Ev + S.Iv; S.Ev = 0; S.Iv = 0;
[~, O] = abm_run(P0, S, 365, C0);
T(end+1,:) = {'No spontaneous infection without infected hosts or vectors', sum(O.G(:,6)) == 0 && O.Ztot == 0, sprintf('%d human, %d mosquito infections', sum(O.G(:,6)), O.Ztot)};
% 5. ITN efficacy: eta_N = 0 increases infection relative to reference
P = P0; P.etaN = 0; set_seed(2); [~, Oa] = abm_run(P, S0, 730, C0); set_seed(2); [~, Ob] = abm_run(P0, S0, 730, C0);
T(end+1,:) = {'ITN protection reduces infection (eta_N = 0 vs 0.70)', sum(Oa.G(:,6)) > sum(Ob.G(:,6)), sprintf('%d vs %d new infections', sum(Oa.G(:,6)), sum(Ob.G(:,6)))};
% 6. treatment: p_tx = 1 shortens symptomatic infection (fewer I_C person-days)
P = P0; P.ptxF = 1; P.ptxM = 1; set_seed(3); [~, Oa] = abm_run(P, S0, 365, C0);
P.ptxF = 0; P.ptxM = 0; set_seed(3); [~, Ob] = abm_run(P, S0, 365, C0);
T(end+1,:) = {'Treatment access lowers prevalence (p_tx = 1 vs 0)', mean(Oa.prev) < mean(Ob.prev), sprintf('%.3f vs %.3f', mean(Oa.prev), mean(Ob.prev))};
% 7. fixed-resource rules allocate exactly floor(0.7N)
nOK = true; for s = {'U','P','D','E','A','H','T'}
    Sx = allocate(P0, S0, struct('alloc', s{1}, 'cov', 0.7, 'Dorder', [1 2 3 4])); nOK = nOK && sum(Sx.z) == floor(0.7*P0.N);
end
T(end+1,:) = {'Each selection rule allocates exactly floor(0.7 N) ITNs', nOK, sprintf('%d ITNs', floor(0.7*P0.N))};
% 8. attribution conserves mosquito infections
[S, O] = abm_run(P0, S0, 365, C0);
T(end+1,:) = {'Attributed events equal new mosquito infections', sum(S.Zhv) + sum(S.ret(:,4)) == O.Ztot, sprintf('%d = %d', sum(S.Zhv) + sum(S.ret(:,4)), O.Ztot)};
fid = fopen('results/A8_verification.txt', 'w');
for i = 1:size(T,1), fprintf(fid, '%-62s %s  (%s)\n', T{i,1}, ifelse(T{i,2}), T{i,3}); end
fclose(fid); type('results/A8_verification.txt');
end
function s = ifelse(b), if b, s = 'PASS'; else, s = 'FAIL'; end, end
