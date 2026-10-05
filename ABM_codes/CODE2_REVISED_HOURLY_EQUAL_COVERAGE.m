%% CODE 2 REVISED: EQUAL NOMINAL COVERAGE VERSUS REALISED PROTECTION
% Standalone stochastic ABM; base MATLAB only.
%
% This revision is aligned with the final hourly Code 1:
%   P0 Equality              : sex-neutral hourly activity, equal use, equal treatment
%   P1 + activity            : literature-constrained hourly activity contrast
%   P2 + effective use       : P1 + female:male effective-use OR = 1.27
%   P3 + treatment barrier   : P2 + counterfactual 10% relative female treatment barrier
%
% IMPORTANT SCIENTIFIC POINTS
% 1. Nominal programme receipt is exactly 70% within each sex in every scenario.
% 2. Human outdoor activity is resolved hourly from 18:00--07:00.
% 3. ITN protection acts on INDOOR exposure only; outdoor exposure is not reduced by ITN use.
% 4. All four scenarios branch from the SAME untreated burn-in state within each replicate.
% 5. The same coverage assignment and matched random-number seed are used across scenarios.
% 6. P3 is a counterfactual stress test, not an empirically estimated sex effect.
% 7. Treatment in P3 is female p_tx = 0.90*baseline; male p_tx remains at baseline.
%
% Figure panels are labelled (a), (b), (c), and (d), as requested.

clear; clc; close all;
tic;

%% --------------------------- CONTROL -----------------------------------
P.N          = 450;
P.R          = 30;
P.burnDays   = 2*365;
P.evalDays   = 5*365;
P.masterSeed = 24092026;

configNames = {'P0 Equality','P1 + activity','P2 + effective use',...
               'P3 + treatment barrier'};
K = numel(configNames);

%% ----------------------- MODEL PARAMETERS ------------------------------
P.dE=12; P.dA=195; P.dC=5; P.dT=5; P.dR=15;
P.pCy=0.60; P.pCa=0.30;
P.qA=0.10; P.qC=1.00; P.qT=0.32;
P.ptx=0.60; P.dSeek=2;

P.a=0.33; P.rhoAge=0.85; P.ageScale=8; P.sigma2e=1.67;
P.m=2.0; P.bvh=0.022; P.bhv=0.48; P.muV=1/14; P.tauV=10.3;
P.xiV=0.20; P.phiV=0;

% Programme
P.etaN     = 0.50;
P.Ceq      = 0.70;
P.piUse0   = 0.85;
P.ORuseFM  = 1.27;
P.deltaRef = -0.15;       % modest male-greater outdoor activity
P.femaleTxMultiplier = 0.90; % P3 counterfactual only

% Hourly exposure structure, identical to final literature-constrained Code 1
P.Rout  = 4.37;
P.hours = [18 19 20 21 22 23 0 1 2 3 4 5 6 7];

P.biteShape = [0.80 1.20 1.10 0.86 0.70 0.92 0.98 0.90 ...
               0.82 0.72 0.62 0.52 0.42 0.32];
P.biteShape = P.biteShape/sum(P.biteShape);

P.pOutOlder = [0.58 0.66 0.63 0.54 0.29 0.15 0.065 0.040 ...
               0.030 0.025 0.025 0.050 0.230 0.500];

P.pOutU5 = [0.40 0.44 0.34 0.20 0.075 0.035 0.015 0.010 ...
            0.008 0.008 0.008 0.018 0.100 0.240];

[pUseF,pUseM] = probabilitiesFromOR(P.piUse0,P.ORuseFM);

fprintf('\n============================================================\n');
fprintf('CODE 2 REVISED: EQUAL COVERAGE AND REALISED PROTECTION\n');
fprintf('============================================================\n');
fprintf('Human agents                    : %d\n',P.N);
fprintf('Matched stochastic runs         : %d\n',P.R);
fprintf('Common untreated burn-in        : %.1f years\n',P.burnDays/365);
fprintf('Evaluation                      : %.1f years\n',P.evalDays/365);
fprintf('Hourly exposure window          : 18:00--07:00 (%d bins)\n',numel(P.hours));
fprintf('Nominal coverage, both sexes    : %.1f%%\n',100*P.Ceq);
fprintf('Reference effective-use mean    : %.1f%%\n',100*P.piUse0);
fprintf('Derived female effective-use p  : %.3f\n',pUseF);
fprintf('Derived male effective-use p    : %.3f\n',pUseM);
fprintf('Female:male effective-use OR    : %.2f\n',P.ORuseFM);
fprintf('Reference activity contrast     : delta_A = %+.2f\n',P.deltaRef);
fprintf('P3 female treatment multiplier  : %.2f (counterfactual)\n',...
        P.femaleTxMultiplier);

%% ---------------------------- STORAGE ----------------------------------
effF=nan(P.R,K); effM=nan(P.R,K);
resF=nan(P.R,K); resM=nan(P.R,K);
recF=nan(P.R,K); recM=nan(P.R,K);
prevF=nan(P.R,K); prevM=nan(P.R,K);
gapClin=nan(P.R,K);

%% ------------------------- MATCHED ANALYSIS -----------------------------
for r=1:P.R
    baseSeed=P.masterSeed+100000*r;

    % One common untreated, sex-neutral burn-in.
    B=burnInCommon(P,baseSeed);

    % One exact 70% coverage assignment per sex, reused in every scenario.
    rng(baseSeed+17,'twister');
    z=false(P.N,1);
    for s=[false true]
        ids=find(B.sex==s);
        nCov=round(P.Ceq*numel(ids));
        z(ids(randperm(numel(ids),nCov)))=true;
    end

    % Same scenario RNG seed across P0--P3 = common random numbers.
    evalSeed=baseSeed+50000;

    for k=1:K
        cfg.deltaA=0;
        cfg.pUseF=P.piUse0; cfg.pUseM=P.piUse0;
        cfg.ptxF=P.ptx; cfg.ptxM=P.ptx;

        if k>=2
            cfg.deltaA=P.deltaRef;
        end
        if k>=3
            cfg.pUseF=pUseF; cfg.pUseM=pUseM;
        end
        if k>=4
            cfg.ptxF=P.ptx*P.femaleTxMultiplier;
            cfg.ptxM=P.ptx; % no artificial male treatment advantage
        end

        O=evaluateCoverageScenario(P,B,z,cfg,evalSeed);

        effF(r,k)=O.effF; effM(r,k)=O.effM;
        resF(r,k)=O.resF; resM(r,k)=O.resM;
        recF(r,k)=O.recF; recM(r,k)=O.recM;
        prevF(r,k)=O.prevF; prevM(r,k)=O.prevM;
        gapClin(r,k)=O.gapClin;
    end
end

%% ----------------------------- REPORT ----------------------------------
fprintf('\n---------------- MAIN MATCHED STOCHASTIC OUTPUT ----------------\n');
fprintf('%-24s EffCovF EffCovM  ResRiskF ResRiskM  RecRiskF RecRiskM  ClinGap/1000py\n',...
    'Configuration');
for k=1:K
    fprintf('%-24s %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f      %8.3f\n',...
        configNames{k},mean(effF(:,k),'omitnan'),mean(effM(:,k),'omitnan'),...
        mean(resF(:,k),'omitnan'),mean(resM(:,k),'omitnan'),...
        mean(recF(:,k),'omitnan'),mean(recM(:,k),'omitnan'),...
        1000*mean(gapClin(:,k),'omitnan'));
end

fprintf('\nFemale-male prevalence differences (percentage points):\n');
for k=1:K
    d=100*(prevF(:,k)-prevM(:,k));
    [lo,hi]=mcCI(d);
    fprintf('%-24s %+7.3f  [%.3f, %.3f]\n',...
        configNames{k},mean(d,'omitnan'),lo,hi);
end

fprintf('\nFemale-male effective-coverage differences (percentage points):\n');
for k=1:K
    d=100*(effF(:,k)-effM(:,k));
    [lo,hi]=mcCI(d);
    fprintf('%-24s %+7.3f  [%.3f, %.3f]\n',...
        configNames{k},mean(d,'omitnan'),lo,hi);
end

%% ----------------------------- FIGURE ----------------------------------
fig=figure('Color','w','Position',[60 60 1300 850]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) realised effective coverage
ax=nexttile; hold(ax,'on');
plotReplicatePairs(ax,effF,effM,configNames,'Realised effective coverage');
ylim(ax,[0.50 0.68]);
title(ax,'(a) Same nominal coverage, realised effective coverage');
grid(ax,'on'); box(ax,'off');

% (b) residual infection risk among programme recipients
ax=nexttile; hold(ax,'on');
plotReplicatePairs(ax,resF,resM,configNames,...
    'P(\geq1 infection | programme recipient)');
ylim(ax,[0 max(0.12,1.15*max([resF(:);resM(:)],[],'omitnan'))]);
title(ax,'(b) Residual infection risk among recipients');
grid(ax,'on'); box(ax,'off');

% (c) recurrent infection risk among programme recipients
ax=nexttile; hold(ax,'on');
plotReplicatePairs(ax,recF,recM,configNames,...
    'P(\geq2 infections | programme recipient)');
ylim(ax,[0 max(0.04,1.20*max([recF(:);recM(:)],[],'omitnan'))]);
title(ax,'(c) Recurrent malaria despite programme receipt');
grid(ax,'on'); box(ax,'off');

% (d) absolute female-male clinical-incidence equity gap
ax=nexttile; hold(ax,'on');
for k=1:K
    y=1000*gapClin(:,k);
    rng(P.masterSeed+900+k,'twister'); % plotting jitter only
    scatter(ax,k+0.12*(rand(P.R,1)-0.5),y,15,'filled',...
        'MarkerFaceAlpha',0.22);
    [lo,hi]=mcCI(y);
    m=mean(y,'omitnan');
    errorbar(ax,k,m,m-lo,hi-m,'k','LineWidth',1.5,'CapSize',8);
end
xticks(ax,1:K); xticklabels(ax,configNames); xtickangle(ax,22);
ylabel(ax,'Absolute female-male clinical-incidence gap per 1000 person-years');
title(ax,'(d) Residual clinical-equity gap');
grid(ax,'on'); box(ax,'off');

title(tl,'Equal nominal coverage does not necessarily imply equal realised protection',...
    'FontWeight','bold');

exportgraphics(fig,'Figure_Code2_EqualCoverage_REVISED.png','Resolution',350);
fprintf('\nFigure saved: Figure_Code2_EqualCoverage_REVISED.png\n');
fprintf('Total runtime: %.1f seconds\n',toc);

%% =======================================================================
% LOCAL FUNCTIONS
% ========================================================================

function B=burnInCommon(P,seed)
rng(seed,'twister');
N=P.N;

sex=false(N,1);
sex(randperm(N,round(N/2)))=true;
age=50*rand(N,1);

e0=exp(sqrt(P.sigma2e)*randn(N,1)-P.sigma2e/2);
ageW=1-P.rhoAge.*exp(-age/P.ageScale);
baseE=e0.*ageW;
baseE=baseE/mean(baseE);

state=zeros(N,1,'uint8');
timer=zeros(N,1,'uint16');
willTreat=false(N,1);
idx=randperm(N,max(2,round(0.05*N)));
state(idx)=2;

Nv0=round(P.m*N);
Sv=round(0.95*Nv0); Ev=round(0.02*Nv0); Iv=Nv0-Sv-Ev;

for t=1:P.burnDays
    pOut=makeHourlyActivityMatrix(P,age,sex,0);
    eEff=hourlyExposureWithIndoorProtection(P,baseE,pOut,false(N,1),0);

    age=age+1/365;
    pClin=P.pCa*ones(N,1); pClin(age<5)=P.pCy;

    sus=state==0;
    newly=sus & rand(N,1)<(1-exp(-P.a*P.bvh*(Iv/N).*eEff));
    state(newly)=1; timer(newly)=0;

    lat=state==1; timer(lat)=timer(lat)+1;
    finish=lat & timer>=P.dE;
    if any(finish)
        rr=rand(N,1);
        clin=finish & rr<pClin; asym=finish & ~clin;
        state(clin)=3; state(asym)=2; timer(finish)=0;
        willTreat(clin)=rand(nnz(clin),1)<P.ptx;
    end

    A=state==2;
    rec=A & rand(N,1)<1-exp(-1/P.dA); state(rec)=5; timer(rec)=0;

    C=state==3; timer(C)=timer(C)+1;
    treat=C & willTreat & timer>=P.dSeek; state(treat)=4; timer(treat)=0;
    C=state==3;
    rec=C & rand(N,1)<1-exp(-1/P.dC);
    state(rec)=5; timer(rec)=0; willTreat(rec)=false;

    T=state==4;
    rec=T & rand(N,1)<1-exp(-1/P.dT); state(rec)=5; timer(rec)=0;

    RR=state==5;
    rec=RR & rand(N,1)<1-exp(-1/P.dR); state(rec)=0; timer(rec)=0;

    c=zeros(N,1);
    c(state==2)=P.qA; c(state==3)=P.qC; c(state==4)=P.qT;
    Ih=sum(eEff.*c)/max(sum(eEff),eps);

    newVE=rbinomFast(Sv,1-exp(-P.a*P.bhv*Ih));
    newVI=rbinomFast(Ev,1-exp(-1/P.tauV));
    dS=rbinomFast(Sv,1-exp(-P.muV));
    dE=rbinomFast(Ev,1-exp(-P.muV));
    dI=rbinomFast(Iv,1-exp(-P.muV));
    recruit=rpoissonFast(max(P.muV*P.m*N*...
        (1+P.xiV*cos(2*pi*(t-P.phiV)/365)),0));

    Sv=max(0,Sv-newVE-dS+recruit);
    Ev=max(0,Ev+newVE-newVI-dE);
    Iv=max(0,Iv+newVI-dI);
end

B.sex=sex; B.age=age; B.baseE=baseE;
B.state=state; B.timer=timer; B.willTreat=willTreat;
B.Sv=Sv; B.Ev=Ev; B.Iv=Iv;
end

function O=evaluateCoverageScenario(P,B,z,cfg,seed)
rng(seed,'twister');
N=P.N;

sex=B.sex; age=B.age; baseE=B.baseE;
state=B.state; timer=B.timer; willTreat=B.willTreat;
Sv=B.Sv; Ev=B.Ev; Iv=B.Iv;

pUse=cfg.pUseM*ones(N,1); pUse(sex)=cfg.pUseF;
pTx=cfg.ptxM*ones(N,1); pTx(sex)=cfg.ptxF;

nInf=zeros(N,1);
useDays=zeros(N,1);
prevFsum=0; prevMsum=0;
clinF=0; clinM=0;

for t=1:P.evalDays
    dayIndex=P.burnDays+t;

    pOut=makeHourlyActivityMatrix(P,age,sex,cfg.deltaA);

    % Daily correct/effective use among nominal programme recipients.
    u=z & (rand(N,1)<pUse);

    % ITN effect applies only to indoor exposure.
    eEff=hourlyExposureWithIndoorProtection(P,baseE,pOut,u,P.etaN);

    age=age+1/365;
    pClin=P.pCa*ones(N,1); pClin(age<5)=P.pCy;

    sus=state==0;
    newly=sus & rand(N,1)<(1-exp(-P.a*P.bvh*(Iv/N).*eEff));
    state(newly)=1; timer(newly)=0;
    nInf(newly)=nInf(newly)+1;
    useDays=useDays+double(u);

    lat=state==1; timer(lat)=timer(lat)+1;
    finish=lat & timer>=P.dE;
    if any(finish)
        rr=rand(N,1);
        clin=finish & rr<pClin; asym=finish & ~clin;
        state(clin)=3; state(asym)=2; timer(finish)=0;
        willTreat(clin)=rand(nnz(clin),1)<pTx(clin);
        clinF=clinF+nnz(clin & sex);
        clinM=clinM+nnz(clin & ~sex);
    end

    A=state==2;
    rec=A & rand(N,1)<1-exp(-1/P.dA); state(rec)=5; timer(rec)=0;

    C=state==3; timer(C)=timer(C)+1;
    treat=C & willTreat & timer>=P.dSeek; state(treat)=4; timer(treat)=0;
    C=state==3;
    rec=C & rand(N,1)<1-exp(-1/P.dC);
    state(rec)=5; timer(rec)=0; willTreat(rec)=false;

    T=state==4;
    rec=T & rand(N,1)<1-exp(-1/P.dT); state(rec)=5; timer(rec)=0;

    RR=state==5;
    rec=RR & rand(N,1)<1-exp(-1/P.dR); state(rec)=0; timer(rec)=0;

    c=zeros(N,1);
    c(state==2)=P.qA; c(state==3)=P.qC; c(state==4)=P.qT;
    Ih=sum(eEff.*c)/max(sum(eEff),eps);

    newVE=rbinomFast(Sv,1-exp(-P.a*P.bhv*Ih));
    newVI=rbinomFast(Ev,1-exp(-1/P.tauV));
    dS=rbinomFast(Sv,1-exp(-P.muV));
    dE=rbinomFast(Ev,1-exp(-P.muV));
    dI=rbinomFast(Iv,1-exp(-P.muV));
    recruit=rpoissonFast(max(P.muV*P.m*N*...
        (1+P.xiV*cos(2*pi*(dayIndex-P.phiV)/365)),0));

    Sv=max(0,Sv-newVE-dS+recruit);
    Ev=max(0,Ev+newVE-newVI-dE);
    Iv=max(0,Iv+newVI-dI);

    infNow=(state==1)|(state==2)|(state==3)|(state==4);
    prevFsum=prevFsum+mean(infNow(sex));
    prevMsum=prevMsum+mean(infNow(~sex));
end

recipF=z & sex; recipM=z & ~sex;

% Realised effective-use coverage:
% population fraction nominally covered AND effectively using protection.
O.effF=P.Ceq*mean(useDays(recipF)/P.evalDays);
O.effM=P.Ceq*mean(useDays(recipM)/P.evalDays);

O.resF=mean(nInf(recipF)>=1);
O.resM=mean(nInf(recipM)>=1);
O.recF=mean(nInf(recipF)>=2);
O.recM=mean(nInf(recipM)>=2);

O.prevF=prevFsum/P.evalDays;
O.prevM=prevMsum/P.evalDays;

pyF=nnz(sex)*P.evalDays/365;
pyM=nnz(~sex)*P.evalDays/365;
JF=clinF/max(pyF,eps);
JM=clinM/max(pyM,eps);
O.gapClin=abs(JF-JM);
end

function pOut=makeHourlyActivityMatrix(P,age,sex,deltaA)
N=numel(age);
base=repmat(P.pOutOlder,N,1);
u5=age<5;
base(u5,:)=repmat(P.pOutU5,nnz(u5),1);

L=logitSafe(base);
shift=zeros(N,1);
shift(sex)=+deltaA/2;
shift(~sex)=-deltaA/2;

pOut=invlogit(L+shift);
pOut=min(max(pOut,0.001),0.95);
end

function eDaily=hourlyExposureWithIndoorProtection(P,baseE,pOut,u,eta)
% Draw hourly outdoor/indoor location.
outdoor=rand(numel(baseE),numel(P.hours))<pOut;

% No-intervention activity exposure used only as a normalization reference.
raw0=(P.Rout*double(outdoor) + double(~outdoor))*P.biteShape(:);

% Protection acts on indoor exposure only.
indoorProtection=1-eta*double(u);
rawP=(P.Rout*double(outdoor) + ...
      double(~outdoor).*indoorProtection)*P.biteShape(:);

% Preserve Code 1's population-scale activity normalization while allowing
% the intervention to lower exposure. Do NOT renormalize rawP itself.
eDaily=baseE.*rawP/max(mean(baseE.*raw0),eps);
end

function plotReplicatePairs(ax,F,M,names,ylab)
K=size(F,2); R=size(F,1);
for k=1:K
    rng(8000+k,'twister'); % jitter only
    xf=k-0.12+0.06*(rand(R,1)-0.5);
    xm=k+0.12+0.06*(rand(R,1)-0.5);

    if k==1
        scatter(ax,xf,F(:,k),13,'filled','MarkerFaceAlpha',0.18,...
            'DisplayName','Female replicates');
        scatter(ax,xm,M(:,k),13,'filled','MarkerFaceAlpha',0.18,...
            'DisplayName','Male replicates');
    else
        scatter(ax,xf,F(:,k),13,'filled','MarkerFaceAlpha',0.18,...
            'HandleVisibility','off');
        scatter(ax,xm,M(:,k),13,'filled','MarkerFaceAlpha',0.18,...
            'HandleVisibility','off');
    end

    [loF,hiF]=mcCI(F(:,k)); mF=mean(F(:,k),'omitnan');
    [loM,hiM]=mcCI(M(:,k)); mM=mean(M(:,k),'omitnan');

    errorbar(ax,k-0.12,mF,mF-loF,hiF-mF,'k','LineWidth',1.3,...
        'CapSize',6,'HandleVisibility','off');
    errorbar(ax,k+0.12,mM,mM-loM,hiM-mM,'k','LineWidth',1.3,...
        'CapSize',6,'HandleVisibility','off');
end
xticks(ax,1:K); xticklabels(ax,names); xtickangle(ax,22);
ylabel(ax,ylab);
legend(ax,'Location','best','Box','off');
end

function [pf,pm]=probabilitiesFromOR(targetMean,OR)
lo=max(1e-5,2*targetMean-1+1e-5);
hi=min(1-1e-5,2*targetMean-1e-5);
for it=1:80
    pm=(lo+hi)/2;
    pf=OR*pm/(1-pm+OR*pm);
    if (pf+pm)/2>targetMean, hi=pm; else, lo=pm; end
end
pm=(lo+hi)/2;
pf=OR*pm/(1-pm+OR*pm);
end

function y=logitSafe(p)
p=min(max(p,1e-8),1-1e-8);
y=log(p./(1-p));
end

function p=invlogit(x)
p=1./(1+exp(-x));
end

function x=rbinomFast(n,p)
if n<=0 || p<=0, x=0; return; end
if p>=1, x=n; return; end
mu=n*p;
if n<=80
    x=sum(rand(n,1)<p);
elseif mu<25
    x=min(n,rpoissonFast(mu));
elseif n*(1-p)<25
    x=n-min(n,rpoissonFast(n*(1-p)));
else
    x=round(mu+sqrt(n*p*(1-p))*randn);
    x=min(max(x,0),n);
end
end

function x=rpoissonFast(lambda)
if lambda<=0, x=0; return; end
if lambda<30
    L=exp(-lambda); k=0; q=1;
    while q>L
        k=k+1; q=q*rand;
    end
    x=k-1;
else
    x=max(0,round(lambda+sqrt(lambda)*randn));
end
end

function [lo,hi]=mcCI(x)
x=x(isfinite(x));
if isempty(x), lo=NaN; hi=NaN; return; end
m=mean(x); h=1.96*std(x)/sqrt(numel(x));
lo=m-h; hi=m+h;
end
