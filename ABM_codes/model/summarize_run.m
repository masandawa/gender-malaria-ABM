function R = summarize_run(P, O)
%SUMMARIZE_RUN  Epidemiological outputs of one simulation (all scalars).
%
%   Prevalence 'prev' counts all blood-stage infections (PCR-like);
%   'rdt' weights non-clinical infections by age-specific detectability.
%   Female:male contrasts are age-standardised over <5 and >=5 years using
%   the simulated population shares.
G  = O.G;  py = G(:,1)/365;                        % person-years, groups fy my fa ma
wA = [sum(G(1:2,1)) sum(G(3:4,1))]/sum(G(:,1));     % age-standardisation weights
stdz = @(v) wA(1)*v([1 2]) + wA(2)*v([3 4]);        % -> [female male]

% ---------------------------------------------------------- population level
R.prev   = sum(G(:,2))/sum(G(:,1));
R.rdt    = sum(G(:,3))/sum(G(:,1));
R.u5rdt  = sum(G(1:2,3))/sum(G(1:2,1));
R.inc    = 1000*sum(G(:,4))/sum(py);               % clinical episodes per 1000 PY
R.incU5  = 1000*sum(G(1:2,4))/sum(py(1:2));
R.nSev   = sum(G(:,5));  R.py = sum(py);            % severe episodes and person-years (pooled later)
R.Ztot   = O.Ztot;                                  % attributed human-to-mosquito infections
R.TI     = double(O.lastActive <= O.nDays - 365);   % no infection in final 365 days

% ---------------------------------------------------------- female-male contrasts
prevG = G(:,2)./G(:,1); rdtG = G(:,3)./G(:,1); J = 1000*G(:,4)./py;
pf = stdz(prevG); rf = stdz(rdtG); jf = stdz(J);
R.prevF = pf(1); R.prevM = pf(2);
R.dPrev = 100*(pf(1) - pf(2));                      % percentage points
R.PRfm  = pf(1)/max(pf(2), eps);
R.PRfmRDT = rf(1)/max(rf(2), eps);
R.PRu5  = prevG(1)/max(prevG(2), eps);              % girls:boys
R.incF  = jf(1); R.incM = jf(2);
R.IRRfm = jf(1)/max(jf(2), eps);
R.Gabs  = wA(1)*abs(J(1) - J(2)) + wA(2)*abs(J(3) - J(4));   % non-cancelling disparity

% ---------------------------------------------------------- adults and pregnancy
SP = O.SP;   % rows: pregnant, non-pregnant women 15-44, men 15-44, all women 15-44
R.PRadult  = (SP(4,2)/SP(4,1))/max(SP(3,2)/SP(3,1), eps);
R.rdtPreg  = SP(1,3)/max(SP(1,1), 1);
R.rdtNonPreg = SP(2,3)/max(SP(2,1), 1);
R.PRpreg   = R.rdtPreg/max(R.rdtNonPreg, eps);       % pregnant : non-pregnant women
R.incPreg  = 1000*SP(1,4)/max(SP(1,1)/365, eps);
R.incNonPreg = 1000*SP(2,4)/max(SP(2,1)/365, eps);
R.IRRpreg  = R.incPreg/max(R.incNonPreg, eps);
R.PregU5   = R.rdtPreg/max(R.u5rdt, eps);           % pregnant : children 6-59 months
R.nSevPreg = SP(1,5); R.pyPreg = SP(1,1)/365;
R.effCovPreg = SP(1,7)/max(SP(1,1), 1);

% ---------------------------------------------------------- exposure and ITN protection
U5 = O.U5; pn = U5(1,2)/max(U5(1,1),1); pu = U5(2,2)/max(U5(2,1),1);
R.ORu5 = (pn/(1-pn))/max(pu/(1-pu), eps);           % odds, children not using vs using ITN
EX = O.EX;
R.outF = EX(1,1)/EX(1,2); R.outM = EX(2,1)/EX(2,2); % outdoor share of exposure
R.piI  = sum(EX(:,3))/sum(EX(:,4));                 % indoor share, ITN non-users
R.piIF = EX(1,3)/EX(1,4); R.piIM = EX(2,3)/EX(2,4);
R.effF = EX(1,5)/EX(1,6); R.effM = EX(2,5)/EX(2,6); % realised effective coverage
R.nomF = EX(1,7)/EX(1,6); R.nomM = EX(2,7)/EX(2,6); % nominal receipt
R.use  = sum(EX(:,5))/sum(EX(:,6));
R.extraNets = O.extraNets;

% ---------------------------------------------------------- cumulative infection risk
% among people present throughout follow-up, by ITN receipt at the start,
% directly standardised over quintiles of expected exposure at the start
R.PresF = NaN; R.PresM = NaN; R.riskRecStd = NaN; R.riskNonStd = NaN; R.RRrecStd = NaN;
if ~isempty(O.z0)
    stay = O.born == O.born0;
    inf1 = O.Ninf >= 1;
    rec  = O.z0 == 1 & stay;
    R.PresF = mean(inf1(rec & O.sex == 0)); R.PresM = mean(inf1(rec & O.sex == 1));
    e = O.eExp0(stay); i1 = inf1(stay); zr = O.z0(stay) == 1;
    q = quint(e); rr = NaN(5,1); rn = NaN(5,1);
    for k = 1:5
        rr(k) = mean(i1(q == k & zr)); rn(k) = mean(i1(q == k & ~zr));
    end
    ok = isfinite(rr) & isfinite(rn);
    R.riskRecStd = mean(rr(ok)); R.riskNonStd = mean(rn(ok));
    R.RRrecStd = R.riskRecStd/max(R.riskNonStd, eps);
end
end

function q = quint(x)
[~, o] = sort(x); n = numel(x); q = zeros(n,1);
q(o) = min(5, 1 + floor(5*(0:n-1)'/n));
end
