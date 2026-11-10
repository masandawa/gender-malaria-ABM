function S = abm_init(P)
%ABM_INIT  Initialise human agents and the aggregate mosquito population.
N = P.N;
S.t = 0;
S.age = trunc_exp(N, P.meanAge, P.maxAge);
S.sex = double(rand(N,1) < 0.5);                 % 0 = female, 1 = male
S.e0  = exp(sqrt(P.sig2e)*randn(N,1) - P.sig2e/2);% persistent propensity (mean 1)
S.re  = P.sigAct*randn(N,1);                     % persistent outdoor-activity effect
% pregnancy (demography is always on; biology is switched by P.pregBio)
S.preg = zeros(N,1); S.pregClock = zeros(N,1); S.iptp = zeros(N,1);
wra = S.sex==0 & S.age>=P.wraMin & S.age<P.wraMax;
fP  = (P.TFR/(P.wraMax-P.wraMin))*P.dPreg/365;  % expected fraction pregnant
idx = find(wra & rand(N,1) < fP);
S.preg(idx) = 1; S.pregClock(idx) = randi(P.dPreg, numel(idx), 1);
% malaria states: 0 S, 1 I_A, 2 I_C, 3 T, 4 R ; latent clock (days left)
S.x   = zeros(N,1);
seed  = rand(N,1) < 0.05; S.x(seed) = 1;
S.lat = zeros(N,1);
S.txFlag = zeros(N,1); S.txTimer = zeros(N,1);
% preventive intervention
S.z = zeros(N,1); S.h = zeros(N,1);           % receipt; habitual use
% individual histories (current occupants of each slot)
S.Ninf = zeros(N,1); S.Nclin = zeros(N,1); S.Nsev = zeros(N,1); S.Zhv = zeros(N,1);
S.winInf = zeros(N,1); S.winZ = zeros(N,1);   % counts since last re-allocation
S.born = zeros(N,1);                           % day the occupant entered
% retired histories (agents replaced through demographic turnover)
S.ret = zeros(0,5);                            % [sex Ninf Nclin Zhv zEver]
S.zEver = zeros(N,1);
% lifetime infection history (used only by the alternative immunity structures);
% initial agents receive histories consistent with their age and exposure
eA = S.e0.*(1 - P.rho*exp(-S.age/P.aage));
S.Nlife = P.lifeInfRate*S.age.*eA;               % expected count (no random draw, so reference streams are unchanged)
S.Dlife = 365*P.lifeInfFrac*S.age.*eA;
% mosquitoes
S.Sv = round(P.m*N); S.Ev = 0; S.Iv = round(0.01*P.m*N);
% activity probabilities (rows updated when age class changes)
S.ageCls = double(S.age < 5);
[S.Pout, S.Bed] = activity_rows(P, S, (1:N)');
end

function a = trunc_exp(n, mu, amax)
a = -mu*log(rand(n,1));
bad = a >= amax;
while any(bad)
    a(bad) = -mu*log(rand(sum(bad),1)); bad = a >= amax;
end
end
