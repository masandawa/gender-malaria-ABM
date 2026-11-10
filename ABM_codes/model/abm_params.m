function P = abm_params()
%ABM_PARAMS  Reference parameter set of the agent-based malaria model.
%   Values, units and sources are listed in Table 2 of the manuscript.
%   Values marked "assumption" are scenario inputs examined in sensitivity
%   analyses. The mosquito-to-human ratio m is calibrated (analyses/A0_calibration.m);
%   the value below is the calibrated value used for the reported results.
%
%   Compatible with MATLAB (R2019b+) and GNU Octave (>= 6).

% ---------------- population and demography -----------------------------
P.N        = 1000;          % human agents
P.meanAge  = 25;            % mean of (truncated) exponential age distribution, years
P.maxAge   = 85;            % truncation, years
P.lifeExp  = 25;            % mean residence/survival -> daily replacement prob
P.TFR      = 4.8;           % total fertility rate, TDHS-MIS 2022
P.wraMin   = 15; P.wraMax = 45;   % reproductive ages, years
P.dPreg    = 280;           % pregnancy duration, days
P.dANC     = 110;           % day of pregnancy at first ANC contact (approx. 16 weeks)

% ---------------- human infection --------------------------------------
P.dE  = 12;   P.dA = 195;  P.dC = 5;  P.dT = 5;  P.dR = 15;     % days
P.pCy = 0.60; P.pCa = 0.30;           % clinical probability <5 and >=5 years
P.qA  = 0.10; P.qC  = 1.00; P.qT = 0.32;   % relative infectiousness
P.ptx = 0.60; P.dseek = 2;            % effective treatment probability, delay (days)
P.ptxF = P.ptx; P.ptxM = P.ptx;       % sex-specific (scenario)
P.pSev = 0.01;                        % probability that a clinical episode is severe (assumption)
% RDT/microscopy detectability of non-clinical infection, declining with age
% because acquired immunity lowers parasite density (Okell et al. 2012):
% theta_A(a) = thetaMin + (thetaMax - thetaMin)*exp(-a/thetaAge)
P.thetaMax = 0.50; P.thetaMin = 0.20; P.thetaAge = 10;

% ---------------- vector and transmission --------------------------------
P.a = 0.33; P.bvh = 0.022; P.bhv = 0.48; P.muv = 1/14; P.tauv = 10.3;
P.m = 19.30;                          % mosquito-to-human ratio (calibrated, A0; overwritten by run_all)
P.xiv = 0.20; P.phiv = 0;             % seasonality

% ---------------- exposure heterogeneity ---------------------------------
P.rho = 0.85; P.aage = 8; P.sig2e = 1.67;
P.sigAct = 0.50;                      % SD of individual random effect in outdoor log-odds (assumption)

% ---------------- hourly activity and biting (18:00-07:00, 14 hours) -----
% relative hourly biting activity (An. arabiensis-dominated profile)
w = [0.04 0.07 0.09 0.10 0.11 0.11 0.10 0.09 0.08 0.07 0.06 0.04 0.03 0.01];
P.w = w / sum(w);
% probability of being outdoors, age >= 5 years and < 5 years
P.poutA = [0.95 0.90 0.80 0.60 0.35 0.15 0.06 0.04 0.04 0.06 0.12 0.30 0.60 0.85];
P.poutY = [0.80 0.60 0.35 0.15 0.05 0.02 0.01 0.01 0.01 0.01 0.03 0.08 0.25 0.50];
% hours spent in bed when indoors (ITN can protect only these hours)
P.bedA  = logical([0 0 0 0 1 1 1 1 1 1 1 0 0 0]);   % 22:00-05:00
P.bedY  = logical([0 0 1 1 1 1 1 1 1 1 1 1 0 0]);   % 20:00-06:00
P.Rout  = 1.00;           % outdoor:indoor biting-rate ratio (mainland Tanzania HLC, ref.)
P.deltaA = -0.15;         % female-male contrast in log-odds of being outdoors (neg. = men more)

% ---------------- preventive intervention (ITN) --------------------------
P.etaN  = 0.70;           % reduction of in-bed exposure while using an ITN
P.piU   = 0.85;           % probability that a recipient is a habitual (persistent) user
P.piUf  = P.piU; P.piUm = P.piU;
P.useNightly = 0;         % 1 = independent nightly use (sensitivity), 0 = habitual use

% ---------------- sex-specific biology (switches) ------------------------
P.pregBio = 0;            % 1 = pregnancy-specific biology on
P.kappaP  = 2.0;          % relative attractiveness of pregnant women (Lindsay et al. 2000)
P.rrClinP = 1.5;          % relative clinical probability during pregnancy (assumption)
P.rrSevP  = 3.0;          % relative risk of severe disease in pregnancy (Desai et al. 2007)
P.iptpOn  = 0;            % 1 = IPTp-SP delivered at ANC
P.cIPTp   = 0.32;         % IPTp3+ coverage, TDHS-MIS 2022
P.effIPTp = 0.40;         % protective efficacy of IPTp-SP against infection (assumption)
P.sexClear = 0;           % 1 = sex-specific clearance of non-clinical infection
P.hrClear  = 1.82;        % female:male clearance hazard ratio (Briggs et al. 2020)
P.clearMinAge = 5;        % sex difference applied from this age (years)

% ---------------- clinical immunity structure (Experiment 6, A9) ----------
% 0 = reference: clinical probability depends on age only (pCy, pCa)
% 1 = infection-count immunity: pC = pCmin + (pCy-pCmin)*exp(-Nlife/immK1)
% 2 = parasite-exposure immunity: pC = pCmin + (pCy-pCmin)*exp(-Dlife/immK2)
%     (Nlife = lifetime infections, Dlife = lifetime days infected)
P.immType = 0; P.pCmin = 0.05; P.immK1 = 18; P.immK2 = 3400;
P.lifeInfRate = 0.55;     % infections per year used only to initialise histories
P.lifeInfFrac = 0.29;     % fraction of time infected used only to initialise histories

% ---------------- calibration target (Morogoro region, TDHS-MIS 2022) ----
P.calTarget = 0.06;       % RDT prevalence, children 6-59 months
P.calLo = 0.04; P.calHi = 0.08;   % acceptance interval (assumed +/-2 percentage points)
end
