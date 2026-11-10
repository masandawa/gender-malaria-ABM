function S = allocate(P, S, C)
%ALLOCATE  Assign ITN receipt z_i under the allocation rule C.alloc.
%   Fixed-resource rules ('U','P','D','E','A','H','T') give exactly floor(C.cov*N)
%   recipients and differ only in the information used to rank individuals.
N = P.N; K = floor(C.cov*N);
zOld = S.z; hOld = S.h;
switch C.alloc
    case 'keep'
        return
    case 'none'
        S.z(:) = 0;
    case 'sex'          % equal nominal coverage within each sex (Analysis 2, P0-P3);
                        % existing recipients keep their ITN, deficits are topped up
        for sx = 0:1
            S.z = top_up(S.z, find(S.sex == sx), C.cov);
        end
    case 'policyProb'   % policy-aligned coverage (TDHS-MIS 2022 informed): pregnant
                        % women, children <5 and others topped up to group targets
        pg = S.preg == 1; u5 = S.age < 5 & ~pg; ot = ~pg & ~u5;
        S.z = top_up(S.z, find(pg), C.covPol(2));
        S.z = top_up(S.z, find(u5), C.covPol(1));
        S.z = top_up(S.z, find(ot), C.covPol(3));
    case 'U'            % uniform (random) selection
        S.z(:) = 0; k = randperm(N); S.z(k(1:K)) = 1;
    case 'P'            % policy-aligned priority: pregnant women, then under-5s, then random
        score = 2*(S.preg == 1) + 1*(S.age < 5 & S.preg == 0);
        S.z = top_k(score, K);
    case 'D'            % demographic: age-sex groups ranked by previous clinical incidence
        g = 1 + S.sex + 2*(S.age >= 5);
        rankG = zeros(4,1); rankG(C.Dorder) = 4:-1:1;
        S.z = top_k(rankG(g), K);
    case 'E'            % exposure-informed: expected activity-adjusted exposure
        om  = 1 - P.rho*exp(-S.age/P.aage);
        kap = ones(N,1); if P.pregBio, kap(S.preg==1) = P.kappaP; end
        Om  = S.Pout*(P.w*P.Rout)' + (1 - S.Pout)*P.w';
        S.z = top_k(S.e0.*om.*kap.*Om, K);
    case 'A'            % observable activity proxy: age, pregnancy and reported night-time
                        % activity, without the unobservable individual attractiveness e0
        om  = 1 - P.rho*exp(-S.age/P.aage);
        kap = ones(N,1); if P.pregBio, kap(S.preg==1) = P.kappaP; end
        Om  = S.Pout*(P.w*P.Rout)' + (1 - S.Pout)*P.w';
        S.z = top_k(om.*kap.*Om, K);
    case 'H'            % infection-history-informed (previous window)
        S.z = top_k(S.winInf, K);
    case 'T'            % transmission-history-informed benchmark (previous window)
        S.z = top_k(S.winZ, K);
    otherwise
        error('Unknown allocation rule %s', C.alloc);
end
S.zEver = double(S.zEver | S.z);
% habitual ITN use is a persistent behaviour drawn when an ITN is received
% (sex-specific); individuals who keep their ITN keep their habit
piUi = P.piUf*(S.sex==0) + P.piUm*(S.sex==1);
hNew = double(rand(N,1) < piUi);
kept = S.z == 1 & zOld == 1;
S.h = zeros(N,1); S.h(kept) = hOld(kept);
nw = S.z == 1 & zOld == 0; S.h(nw) = hNew(nw);
S.winInf(:) = 0; S.winZ(:) = 0;
end

function z = top_k(score, K)
% select the K highest scores; ties broken at random
n = numel(score);
[~, o] = sortrows([-score(:) rand(n,1)]);
z = zeros(n,1); z(o(1:K)) = 1;
end

function z = top_up(z, idx, cov)
% keep existing recipients in idx, then add (or remove) at random to reach cov
target = round(cov*numel(idx));
has = idx(z(idx) == 1); hasnt = idx(z(idx) == 0);
if numel(has) > target
    has = has(randperm(numel(has))); z(has(target+1:end)) = 0;
elseif numel(has) < target
    hasnt = hasnt(randperm(numel(hasnt))); z(hasnt(1:target-numel(has))) = 1;
end
end
