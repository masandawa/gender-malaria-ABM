function [S, O] = abm_run(P, S, nDays, C)
%ABM_RUN  Simulate the agent-based model for nDays daily steps.
%
%   [S, O] = abm_run(P, S, nDays, C) advances state S under parameters P and
%   control structure C, returning the final state and output accumulators O.
%
%   C fields (all optional):
%     alloc    'none' | 'keep' | 'sex' | 'policyProb' | 'U' | 'P' | 'D' | 'E' | 'A' | 'H' | 'T'
%     cov      nominal coverage (default 0.70)
%     covPol   [under-5, pregnant, other] receipt probabilities for 'policyProb'
%     ancTopUp probability that a pregnant woman without an ITN receives one at first ANC
%     realloc  re-allocation interval in days (Inf = allocate once at day 1)
%     Dorder   age-sex group priority order for demographic allocation
%     resetHist  1 = reset individual histories at day 1 (evaluation period)
%     redrawHabit 1 = draw ITN-use habits of recipients from P at day 1
%
%   Group index g: 1 female <5, 2 male <5, 3 female >=5, 4 male >=5.

if nargin < 4, C = struct(); end
C = set_default(C, 'alloc', 'keep');   C = set_default(C, 'cov', 0.70);
C = set_default(C, 'covPol', [0.78 0.78 0.68]);
C = set_default(C, 'ancTopUp', 0);     C = set_default(C, 'realloc', Inf);
C = set_default(C, 'Dorder', [1 2 3 4]); C = set_default(C, 'resetHist', 1);
C = set_default(C, 'redrawHabit', 1);

N  = P.N;
PA = 1-exp(-1/P.dA); PC = 1-exp(-1/P.dC); PT = 1-exp(-1/P.dT); PR = 1-exp(-1/P.dR);
PEI = 1-exp(-1/P.tauv); PD = 1-exp(-P.muv);
pDie = 1/(P.lifeExp*365);
pConc = (P.TFR/(P.wraMax-P.wraMin))/365 / (1 - (P.TFR/(P.wraMax-P.wraMin))*P.dPreg/365);
wR = P.w*P.Rout;
Lam0 = P.muv*P.m*N;

% ------------------------------------------------------------------ outputs
O = init_outputs(nDays);
if C.resetHist
    S.Ninf(:) = 0; S.Nclin(:) = 0; S.Nsev(:) = 0; S.Zhv(:) = 0; S.ret = zeros(0,5);
    S.zEver(:) = 0;
end
tStart = S.t;

for d = 1:nDays
    S.t = S.t + 1; t = S.t;
    % ---------------------------------------------------------- allocation
    if d == 1 || (isfinite(C.realloc) && mod(d-1, C.realloc) == 0)
        S = allocate(P, S, C);
        if d == 1
            if C.redrawHabit && ~strcmp(C.alloc, 'keep')
                % habits of current recipients follow the scenario's use parameters
                piUi = P.piUf*(S.sex==0) + P.piUm*(S.sex==1);
                S.h = double(S.z == 1 & rand(N,1) < piUi);
            end
            O.z0 = S.z; O.Ninf0 = S.Ninf; O.born0 = S.born;
            om0 = 1 - P.rho*exp(-S.age/P.aage);
            O.eExp0 = S.e0.*om0.*(S.Pout*(P.w*P.Rout)' + (1 - S.Pout)*P.w');
        end
    end
    % --------------------------------------------------------- demography
    S.age = S.age + 1/365;
    die = find(rand(N,1) < pDie | S.age >= P.maxAge);
    if ~isempty(die)
        S.ret = [S.ret; S.sex(die) S.Ninf(die) S.Nclin(die) S.Zhv(die) S.zEver(die)];
        n = numel(die);
        S.age(die) = 0; S.sex(die) = double(rand(n,1) < 0.5);
        S.e0(die)  = exp(sqrt(P.sig2e)*randn(n,1) - P.sig2e/2);
        S.re(die)  = P.sigAct*randn(n,1);
        S.preg(die) = 0; S.pregClock(die) = 0; S.iptp(die) = 0;
        S.x(die) = 0; S.lat(die) = 0; S.txFlag(die) = 0; S.txTimer(die) = 0;
        S.z(die) = 0; S.zEver(die) = 0; S.h(die) = 0;
        S.Ninf(die) = 0; S.Nclin(die) = 0; S.Nsev(die) = 0; S.Zhv(die) = 0;
        S.winInf(die) = 0; S.winZ(die) = 0; S.born(die) = t;
        S.Nlife(die) = 0; S.Dlife(die) = 0;
    end
    newCls = double(S.age < 5);
    dm = false(N,1); dm(die) = true;
    chg = find(newCls ~= S.ageCls | dm);
    if ~isempty(chg)
        S.ageCls(chg) = newCls(chg);
        [S.Pout(chg,:), S.Bed(chg,:)] = activity_rows(P, S, chg);
    end
    % ---------------------------------------------------------- pregnancy
    pg = S.preg == 1;
    S.pregClock(pg) = S.pregClock(pg) + 1;
    fin = pg & S.pregClock >= P.dPreg;
    S.preg(fin) = 0; S.pregClock(fin) = 0; S.iptp(fin) = 0;
    elig = S.sex == 0 & S.preg == 0 & S.age >= P.wraMin & S.age < P.wraMax;
    conc = elig & rand(N,1) < pConc;
    S.preg(conc) = 1; S.pregClock(conc) = 0;
    if P.iptpOn, S.iptp(conc) = double(rand(sum(conc),1) < P.cIPTp); end
    if C.ancTopUp > 0
        anc = S.preg == 1 & S.pregClock == P.dANC & S.z == 0 & rand(N,1) < C.ancTopUp;
        S.z(anc) = 1; S.zEver(anc) = 1; O.extraNets = O.extraNets + sum(anc);
        S.h(anc) = double(rand(sum(anc),1) < P.piUf);
    end
    % ------------------------------------------------- activity and exposure
    Lout = rand(N,14) < S.Pout;
    inBed = ~Lout & S.Bed;  inAw = ~Lout & ~S.Bed;
    om  = 1 - P.rho*exp(-S.age/P.aage);
    kap = ones(N,1);
    if P.pregBio, kap(S.preg==1) = P.kappaP; end
    eT  = S.e0 .* om .* kap;
    eOut = eT .* (double(Lout)*wR');
    eAw  = eT .* (double(inAw)*P.w');
    eBed = eT .* (double(inBed)*P.w');
    if P.useNightly                              % sensitivity: independent nightly use
        piUi = P.piUf*(S.sex==0) + P.piUm*(S.sex==1);
        use  = S.z == 1 & rand(N,1) < piUi;
    else                                         % reference: habitual use drawn at receipt
        use  = S.z == 1 & S.h == 1;
    end
    eAct = eOut + eAw + eBed;
    eEff = eOut + eAw + eBed .* (1 - P.etaN*use);
    % ------------------------------------------------ mosquito-to-human
    x0 = S.x; lat0 = S.lat;
    lam = P.a*P.bvh*S.Iv/N .* eEff;
    if P.iptpOn, lam = lam .* (1 - P.effIPTp*(S.iptp==1 & S.preg==1)); end
    susc = x0 == 0 & lat0 == 0;
    newInf = susc & rand(N,1) < (1-exp(-lam));
    % ------------------------------------------------ within-host progression
    S.lat(lat0 > 0) = lat0(lat0 > 0) - 1;
    done = find(lat0 == 1);
    clinNew = false(N,1); sevNew = false(N,1);
    if ~isempty(done)
        if P.immType == 1
            pC = P.pCmin + (P.pCy - P.pCmin)*exp(-S.Nlife(done)/P.immK1);
        elseif P.immType == 2
            pC = P.pCmin + (P.pCy - P.pCmin)*exp(-S.Dlife(done)/P.immK2);
        else
            pC = P.pCa*ones(numel(done),1); pC(S.age(done) < 5) = P.pCy;
        end
        ad5 = S.age(done) >= 5;
        O.PC = O.PC + [sum(pC(ad5 & S.sex(done)==0)) sum(ad5 & S.sex(done)==0); ...
                       sum(pC(ad5 & S.sex(done)==1)) sum(ad5 & S.sex(done)==1)];
        isP = S.preg(done) == 1 & P.pregBio;
        pC(isP) = min(1, pC(isP)*P.rrClinP);
        cl = rand(numel(done),1) < pC;
        dc = done(cl); da = done(~cl);
        S.x(da) = 1;
        S.x(dc) = 2; S.Nclin(dc) = S.Nclin(dc) + 1; clinNew(dc) = true;
        ptxI = P.ptxF*(S.sex(dc)==0) + P.ptxM*(S.sex(dc)==1);
        S.txFlag(dc) = double(rand(numel(dc),1) < ptxI); S.txTimer(dc) = P.dseek;
        pS = P.pSev*ones(numel(dc),1);
        if P.pregBio, pS(S.preg(dc)==1) = min(1, P.pSev*P.rrSevP); end
        sv = dc(rand(numel(dc),1) < pS); S.Nsev(sv) = S.Nsev(sv) + 1; sevNew(sv) = true;
    end
    % clinical (existing at start of day)
    ic = find(x0 == 2);
    if ~isempty(ic)
        S.txTimer(ic) = S.txTimer(ic) - 1;
        tr = ic(S.txFlag(ic) == 1 & S.txTimer(ic) <= 0);
        S.x(tr) = 3;
        un = ic(S.txFlag(ic) == 0);
        S.x(un(rand(numel(un),1) < PC)) = 1;    % untreated: symptoms resolve, parasites persist
    end
    % non-clinical clearance (optionally sex-specific)
    ia = find(x0 == 1);
    if ~isempty(ia)
        pa = PA*ones(numel(ia),1);
        if P.sexClear
            ad = S.age(ia) >= P.clearMinAge;
            dF = P.dA/sqrt(P.hrClear); dM = P.dA*sqrt(P.hrClear);
            pa(ad & S.sex(ia)==0) = 1-exp(-1/dF);
            pa(ad & S.sex(ia)==1) = 1-exp(-1/dM);
        end
        clr = ia(rand(numel(ia),1) < pa);
        S.x(clr) = 4;
        a5 = S.age(clr) >= 5;
        O.IAX = O.IAX + [sum(a5 & S.sex(clr)==0) sum(a5 & S.sex(clr)==1)];
    end
    it = find(x0 == 3); S.x(it(rand(numel(it),1) < PT)) = 4;
    ir = find(x0 == 4); S.x(ir(rand(numel(ir),1) < PR)) = 0;
    ni = find(newInf);
    S.lat(ni) = P.dE; S.Ninf(ni) = S.Ninf(ni) + 1; S.winInf(ni) = S.winInf(ni) + 1;
    S.Nlife(ni) = S.Nlife(ni) + 1;
    S.Dlife = S.Dlife + double(x0 >= 1 & x0 <= 3);
    % ------------------------------------------------ human-to-mosquito
    c = P.qA*(x0==1) + P.qC*(x0==2) + P.qT*(x0==3);
    v = eEff .* c; sv = sum(v); se = sum(eEff);
    BSE = 0;
    if sv > 0 && S.Sv > 0
        Ih  = sv/se;
        BSE = rbinom(S.Sv, 1-exp(-P.a*P.bhv*Ih));
        if BSE > 0
            cs = cumsum(v);
            r  = rand(BSE,1)*cs(end);
            id = lookup_idx(cs, r);
            cnt = accumarray(id, 1, [N 1]);
            S.Zhv = S.Zhv + cnt; S.winZ = S.winZ + cnt;
        end
    end
    % ------------------------------------------------ mosquito dynamics
    BEI = rbinom(S.Ev, PEI);
    S.Sv = S.Sv - BSE; S.Ev = S.Ev + BSE - BEI; S.Iv = S.Iv + BEI;
    S.Sv = S.Sv - rbinom(S.Sv, PD); S.Ev = S.Ev - rbinom(S.Ev, PD); S.Iv = S.Iv - rbinom(S.Iv, PD);
    S.Sv = S.Sv + rpois(Lam0*(1 + P.xiv*cos(2*pi*(t - P.phiv)/365)));
    % ------------------------------------------------ recording
    inf  = double(x0 >= 1 & x0 <= 3);
    thA  = P.thetaMin + (P.thetaMax - P.thetaMin)*exp(-S.age/P.thetaAge);
    rdt  = double(x0==2 | x0==3) + thA.*double(x0==1);
    yng  = S.age < 5; fem = S.sex == 0;
    Gm   = double([fem & yng, ~fem & yng, fem & ~yng, ~fem & ~yng]);
    M    = [ones(N,1) inf rdt double(clinNew) double(sevNew) double(newInf) double(use)];
    O.G  = O.G + Gm'*M;                          % 4 groups x [pd inf rdt clin sev ninf use]
    wra  = fem & S.age>=15 & S.age<45;
    sp   = double([S.preg==1, wra & S.preg==0, ~fem & S.age>=15 & S.age<45, wra]);
    O.SP = O.SP + sp'*M;    % pregnant, non-pregnant women 15-44, men 15-44, all women 15-44
    u5 = S.age >= 0.5 & yng;
    O.U5 = O.U5 + [sum(u5 & ~use) sum(rdt(u5 & ~use)); sum(u5 & use) sum(rdt(u5 & use))];
    nu = double(~use);
    Sx = double([fem, ~fem]);
    uu = double(use);
    O.EX = O.EX + Sx'*[eOut eAct (eAw+eBed).*nu eAct.*nu uu ones(N,1) S.z (eAw+eBed.*(1-P.etaN)).*uu eEff.*uu];
    O.Ztot = O.Ztot + BSE;
    O.prev(d) = mean(inf); O.rdtAll(d) = mean(rdt);
    O.SEV = O.SEV + [sum(sevNew & fem) sum(sevNew & ~fem)];
    st = double([x0==0 & lat0==0, lat0>0, x0==1, x0==2, x0==3, x0==4]);
    ad = ~yng;
    O.ST = O.ST + [double(fem & ad)'; double(~fem & ad)']*st;   % adults >=5 y: person-days by state
    if any(inf) || any(S.lat > 0) || (S.Ev + S.Iv) > 0, O.lastActive = d; end
end
O.nDays = nDays; O.tStart = tStart;
O.Ninf = S.Ninf; O.Nclin = S.Nclin; O.Zhv = S.Zhv; O.sex = S.sex; O.ret = S.ret;
O.z = S.z; O.born = S.born; O.age = S.age;
end

% =========================================================================
function O = init_outputs(nDays)
% O.G : 4 age-sex groups x [person-days infected rdt clinical severe newInf ITN-use]
% O.SP: pregnant women, non-pregnant women 15-44, men 15-44, all women 15-44 x same columns
% O.EX: sex (F,M) x [eOut eAct eIndoor(non-users) eAct(non-users) use-days person-days receipt-days
%       eIndoor(users, residual) eEff(users)]
O.G = zeros(4,7); O.SP = zeros(4,7); O.U5 = zeros(2,2); O.EX = zeros(2,9);
O.SEV = [0 0]; O.Ztot = 0; O.prev = zeros(nDays,1); O.rdtAll = zeros(nDays,1); O.lastActive = 0;
O.extraNets = 0; O.z0 = []; O.Ninf0 = []; O.born0 = []; O.eExp0 = [];
O.ST = zeros(2,6); O.IAX = [0 0]; O.PC = zeros(2,2);   % checks of the clearance pathway (A9)
end

function C = set_default(C, f, v)
if ~isfield(C, f), C.(f) = v; end
end

function id = lookup_idx(cs, r)
% index of the first cumulative weight >= r (inverse-CDF sampling);
% individuals with zero weight can never be selected.
id = sum(bsxfun(@lt, cs(:), r(:)'), 1)' + 1;
id = max(1, min(numel(cs), id));
end
