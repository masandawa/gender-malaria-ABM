%% CODE 3: FIXED-RESOURCE TARGETING AND STOCHASTIC TRANSMISSION INTERRUPTION
% Standalone stochastic ABM analysis
% Requires base MATLAB only.
%
% Five allocation rules, all constrained to the same C_R resource envelope:
%   U = uniform random
%   D = age-sex demographic burden targeting
%   E = individual exposure-informed targeting
%   H = observed infection-history targeting
%   T = individual transmission-history targeting
%
% Outputs:
%   - mean malaria prevalence
%   - clinical incidence
%   - female-male clinical equity gap
%   - cumulative human-to-mosquito infections
%   - sustained local stochastic transmission-interruption probability
%
% Transmission-history targeting uses perfectly observed individual H-to-M
% attribution during the pre-policy history window and is therefore treated
% as an information-rich benchmark rather than routine surveillance data.

clear; clc; close all;
tic;

%% --------------------------- CONTROL -----------------------------------
P.N          = 750;
P.R          = 100;
P.burnDays   = 3*365;
P.historyDays= 365;
P.followDays = 5*365;
P.T0         = 365;
P.masterSeed = 27092026;

P.CR         = 0.70;    % identical resource envelope
P.piUse      = 0.85;
P.etaN       = 0.50;

strategyNames={'Uniform','Demographic','Exposure','Clinical history','Transmission history'};
K=numel(strategyNames);

%% ------------------------- MODEL VALUES --------------------------------
P.dE=12; P.dA=195; P.dC=5; P.dT=5; P.dR=15;
P.pCy=.60; P.pCa=.30;
P.qA=.10; P.qC=1.00; P.qT=.32;
P.ptx=.60; P.dSeek=2;

P.a=.33; P.rhoAge=.85; P.ageScale=8; P.sigma2e=1.67;
P.m=2; P.bvh=.022; P.bhv=.48; P.muV=1/14; P.tauV=10.3;
P.xiV=.20; P.phiV=0;

P.Rout=4.37;
P.deltaRef=-0.15;
P.hours=[18 19 20 21 22 23 0 1 2 3 4 5 6 7];
P.biteShape=[0.80 1.20 1.10 0.86 0.70 0.92 0.98 0.90 0.82 0.72 0.62 0.52 0.42 0.32];
P.biteShape=P.biteShape/sum(P.biteShape);
P.pOutOlder=[0.58 0.66 0.63 0.54 0.29 0.15 0.065 0.040 0.030 0.025 0.025 0.050 0.230 0.500];
P.pOutU5=[0.40 0.44 0.34 0.20 0.075 0.035 0.015 0.010 0.008 0.008 0.008 0.018 0.100 0.240];

fprintf('\n============================================================\n');
fprintf('CODE 3 REVISED: FIXED-RESOURCE TARGETING + STOCHASTIC INTERRUPTION\n');
fprintf('============================================================\n');
fprintf('N humans                  : %d\n',P.N);
fprintf('Monte Carlo realizations  : %d\n',P.R);
fprintf('Resource envelope         : %.0f%%\n',100*P.CR);
fprintf('History window            : %d days\n',P.historyDays);
fprintf('Follow-up                 : %.1f years\n',P.followDays/365);
fprintf('Interruption window T0      : %d days\n',P.T0);

%% -------------------------- STORAGE ------------------------------------
prev=nan(P.R,K);
clinInc=nan(P.R,K);
eqGap=nan(P.R,K);
trans=nan(P.R,K);
extinct=false(P.R,K);
tExt=nan(P.R,K);

%% ----------------------- MATCHED REPLICATES ----------------------------
for r=1:P.R
    % One common pre-policy population/history per replicate
    baseSeed=P.masterSeed+r;
    B=generateBaseline(P,baseSeed);

    % allocation masks created from exactly the same baseline agents
    Z=cell(1,K);
    Z{1}=selectUniform(P,B,baseSeed+100000);
    Z{2}=selectDemographic(P,B);
    Z{3}=selectTop(B.exposureScore,P.CR,baseSeed+200001);
    Z{4}=selectTop(B.nInfHistory,P.CR,baseSeed+200002);
    Z{5}=selectTop(B.ZhvHistory,P.CR,baseSeed+200003);

    for k=1:K
        % same post-policy master seed across strategies within replicate
        O=runPostPolicy(P,B,Z{k},baseSeed+500000);
        prev(r,k)=O.prev;
        clinInc(r,k)=O.clinInc;
        eqGap(r,k)=O.eqGap;
        trans(r,k)=O.trans;
        extinct(r,k)=O.extinct;
        tExt(r,k)=O.tExt;
    end
end

%% ---------------------------- REPORT -----------------------------------
fprintf('\n---------------- POLICY OUTPUT ----------------\n');
fprintf('%-24s Prev%%   Clin/1000py  EqGap/1000py  H->M infections  Interruption%%\n',...
    'Strategy');
for k=1:K
    fprintf('%-24s %6.2f     %8.2f       %8.2f       %10.1f      %7.1f\n',...
        strategyNames{k},100*mean(prev(:,k),'omitnan'),...
        1000*mean(clinInc(:,k),'omitnan'),...
        1000*mean(eqGap(:,k),'omitnan'),...
        mean(trans(:,k),'omitnan'),100*mean(extinct(:,k)));
end

% Use ONE stable uniform denominator.  Do not divide each paired difference
% by its own stochastic uniform count, because near-zero uniform counts can
% create extreme and misleading percentage ratios.
meanTransU=mean(trans(:,1),'omitnan');

fprintf('\nPaired change from uniform allocation:\n');
for k=2:K
    dp=100*(prev(:,k)-prev(:,1));
    dtr=100*(trans(:,1)-trans(:,k))/max(meanTransU,eps);
    fprintf('%-24s prevalence %+7.3f pp; H->M reduction %+7.2f%%; interruption %+6.1f pp\n',...
        strategyNames{k},mean(dp,'omitnan'),mean(dtr,'omitnan'),...
        100*(mean(extinct(:,k))-mean(extinct(:,1))));
end

%% ---------------------------- FIGURE -----------------------------------
fig=figure('Color','w','Position',[45 45 1450 900]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) prevalence
ax=nexttile; hold(ax,'on');
dotWhisker(ax,prev*100,strategyNames,'Malaria prevalence (%)');
title(ax,'(a) Post-allocation malaria prevalence','FontWeight','bold');
formatPanel(ax);

% (b) clinical incidence
ax=nexttile; hold(ax,'on');
dotWhisker(ax,clinInc*1000,strategyNames,...
    'Clinical episodes per 1000 person-years');
title(ax,'(b) Clinical burden','FontWeight','bold');
formatPanel(ax);

% (c) paired H->M transmission reduction relative to uniform
% The common denominator prevents unstable percentages when a particular
% uniform replicate has very few transmission events.
ax=nexttile; hold(ax,'on');
red=zeros(P.R,K);
for k=2:K
    red(:,k)=100*(trans(:,1)-trans(:,k))/max(meanTransU,eps);
end
dotWhisker(ax,red,strategyNames,'H-to-M transmission reduction vs uniform (%)');
yline(ax,0,'--','Color',[0.35 0.35 0.35],'LineWidth',1.1);
title(ax,'(c) Human-to-mosquito transmission','FontWeight','bold');
formatPanel(ax);

% Keep the visible range driven by the central stochastic distribution,
% while preserving all mean/CI summaries.
allRed=red(:,2:end); allRed=allRed(isfinite(allRed));
if ~isempty(allRed)
    qlo=localPercentile(allRed,2.5);
    qhi=localPercentile(allRed,97.5);
    pad=max(5,0.15*(qhi-qlo));
    ylim(ax,[min(-5,qlo-pad) max(10,qhi+pad)]);
end

% (d) sustained local stochastic transmission interruption
ax=nexttile; hold(ax,'on');
pExt=mean(extinct,1);
se=sqrt(pExt.*(1-pExt)/P.R);
b=bar(ax,1:K,100*pExt,0.62,'FaceColor',[0.30 0.55 0.78],...
    'EdgeColor',[0.15 0.15 0.15],'LineWidth',0.7);
errorbar(ax,1:K,100*pExt,100*1.96*se,'k.',...
    'LineWidth',1.5,'CapSize',9);
xticks(ax,1:K); xticklabels(ax,strategyNames); xtickangle(ax,18);
ylabel(ax,'Probability of sustained interruption (%)');
ylim(ax,[0 100]);
title(ax,sprintf('(d) Infection-free for \\geq %d consecutive days',P.T0),...
    'FontWeight','bold');
formatPanel(ax);


exportgraphics(fig,'Figure_Code3_Targeting_PUBLICATION.png','Resolution',400);
fprintf('\nFigure saved: Figure_Code3_Targeting_PUBLICATION.png\n');
fprintf('Runtime: %.1f seconds\n',toc);

%% =======================================================================
% LOCAL FUNCTIONS
% ========================================================================

function B=generateBaseline(P,seed)
rng(seed,'twister');
N=P.N;

sex=false(N,1); sex(randperm(N,round(N/2)))=true;
age=50*rand(N,1);
e0=exp(sqrt(P.sigma2e)*randn(N,1)-P.sigma2e/2);
baseE=e0.*(1-P.rhoAge.*exp(-age/P.ageScale));
baseE=baseE/mean(baseE);


state=zeros(N,1,'uint8'); timer=zeros(N,1,'uint16'); willTreat=false(N,1);
idx0=randperm(N,max(2,round(.05*N))); state(idx0)=2;
Nv0=round(P.m*N); Sv=round(.95*Nv0); Ev=round(.02*Nv0); Iv=Nv0-Sv-Ev;

nInfHist=zeros(N,1); nClinHist=zeros(N,1); Zhist=zeros(N,1);
expoAccum=zeros(N,1); expoDays=0;

for t=1:(P.burnDays+P.historyDays)
    [state,timer,willTreat,Sv,Ev,Iv,events,eNow,age]=...
        oneDay(P,state,timer,willTreat,Sv,Ev,Iv,age,baseE,sex,...
               false(N,1),P.piUse,ones(N,1)*P.ptx,t);
    if t>P.burnDays
        nInfHist=nInfHist+events.newInf;
        nClinHist=nClinHist+events.newClin;
        Zhist=Zhist+events.Zhv;
        expoAccum=expoAccum+eNow;
        expoDays=expoDays+1;
    end
end

% group burden score for demographic strategy
grp=groupIndex(sex,age);
grpBurden=zeros(4,1);
for g=1:4
    grpBurden(g)=sum(nClinHist(grp==g));
end

B.sex=sex; B.age=age; B.baseE=baseE;
B.state=state; B.timer=timer; B.willTreat=willTreat;
B.Sv=Sv; B.Ev=Ev; B.Iv=Iv;
B.nInfHistory=nInfHist;
B.nClinHistory=nClinHist;
B.ZhvHistory=Zhist;
B.exposureScore=expoAccum/max(expoDays,1);
B.group=grp;
B.groupBurden=grpBurden;
end

function z=selectUniform(P,B,seed)
rng(seed,'twister');
N=numel(B.sex); n=round(P.CR*N);
z=false(N,1); z(randperm(N,n))=true;
end

function z=selectDemographic(P,B)
N=numel(B.sex); n=round(P.CR*N);
score=B.groupBurden(B.group);
% add deterministic tiny individual exposure tie-breaker
score=score+1e-8*B.exposureScore;
[~,ord]=sort(score,'descend');
z=false(N,1); z(ord(1:n))=true;
end

function z=selectTop(score,C,seed)
N=numel(score); n=round(C*N);
% Random tie-breaking avoids index-order bias when agents have identical scores.
rng(seed,'twister');
score=double(score(:));
jitter=1e-10*randn(N,1);
[~,ord]=sort(score+jitter,'descend');
z=false(N,1); z(ord(1:n))=true;
end

function O=runPostPolicy(P,B,z,seed)
rng(seed,'twister');
N=numel(B.sex);

sex=B.sex; age=B.age; baseE=B.baseE;
state=B.state; timer=B.timer; willTreat=B.willTreat;
Sv=B.Sv; Ev=B.Ev; Iv=B.Iv;

prevSum=0; nPrevDays=0;
nClin=0; clinF=0; clinM=0;
pyF=0; pyM=0;
transTotal=0;

zeroRun=0; firstExt=nan;

for t=1:P.followDays
    [state,timer,willTreat,Sv,Ev,Iv,events,~,age]=...
        oneDay(P,state,timer,willTreat,Sv,Ev,Iv,age,baseE,sex,...
               z,P.piUse,ones(N,1)*P.ptx,t);

    infNow=state==1|state==2|state==3|state==4;
    prevSum=prevSum+mean(infNow);
    nPrevDays=nPrevDays+1;

    nClin=nClin+sum(events.newClin);
    clinF=clinF+sum(events.newClin & sex);
    clinM=clinM+sum(events.newClin & ~sex);
    pyF=pyF+nnz(sex)/365; pyM=pyM+nnz(~sex)/365;
    transTotal=transTotal+sum(events.Zhv);

    if ~any(infNow) && Ev==0 && Iv==0
        zeroRun=zeroRun+1;
        if zeroRun>=P.T0 && isnan(firstExt)
            firstExt=t-P.T0+1;
        end
    else
        zeroRun=0;
    end
end

O.prev=prevSum/max(nPrevDays,1);
O.clinInc=nClin/(N*P.followDays/365);
JF=clinF/max(pyF,eps); JM=clinM/max(pyM,eps);
O.eqGap=abs(JF-JM);
O.trans=transTotal;
O.extinct=~isnan(firstExt);
O.tExt=firstExt;
end

function [state,timer,willTreat,Sv,Ev,Iv,E,eEff,age]=...
    oneDay(P,state,timer,willTreat,Sv,Ev,Iv,age,baseE,sex,z,pUse,pTx,t)

N=numel(state);
age=age+1/365;
pClin=P.pCa*ones(N,1); pClin(age<5)=P.pCy;

% Literature-constrained hourly activity, aligned with final Code 1.
pOut=makeHourlyActivityMatrix(P,age,sex,P.deltaRef);
outdoor=rand(N,numel(P.hours))<pOut;

% Baseline hourly exposure and indoor-only programme protection.
raw0=(P.Rout*double(outdoor)+double(~outdoor))*P.biteShape(:);
u=z & (rand(N,1)<pUse);
indoorProtection=1-P.etaN*double(u);
rawP=(P.Rout*double(outdoor)+double(~outdoor).*indoorProtection)*P.biteShape(:);

% Normalize to the no-programme activity scale, not the protected scale.
% This preserves genuine intervention-induced reductions in exposure.
eEff=baseE.*rawP/max(mean(baseE.*raw0),eps);

E.newInf=zeros(N,1); E.newClin=zeros(N,1); E.Zhv=zeros(N,1);

sus=state==0;
pH=1-exp(-P.a*P.bvh*(Iv/N).*eEff);
newly=sus & (rand(N,1)<pH);
state(newly)=1; timer(newly)=0; E.newInf=double(newly);

lat=state==1; timer(lat)=timer(lat)+1;
finish=lat & timer>=P.dE;
if any(finish)
    rr=rand(N,1); clin=finish & rr<pClin; asym=finish & ~clin;
    state(clin)=3; state(asym)=2; timer(finish)=0;
    willTreat(clin)=rand(nnz(clin),1)<pTx(clin);
    E.newClin=double(clin);
end

A=state==2;
rec=A & rand(N,1)<1-exp(-1/P.dA); state(rec)=5; timer(rec)=0;

C=state==3; timer(C)=timer(C)+1;
treat=C & willTreat & timer>=P.dSeek; state(treat)=4; timer(treat)=0;
C=state==3;
rec=C & rand(N,1)<1-exp(-1/P.dC);
state(rec)=5; timer(rec)=0; willTreat(rec)=false;

TT=state==4;
rec=TT & rand(N,1)<1-exp(-1/P.dT); state(rec)=5; timer(rec)=0;

RR=state==5;
rec=RR & rand(N,1)<1-exp(-1/P.dR); state(rec)=0; timer(rec)=0;

c=zeros(N,1);
c(state==2)=P.qA; c(state==3)=P.qC; c(state==4)=P.qT;
w=eEff.*c;
Ih=sum(w)/max(sum(eEff),eps);
newVE=rbinomFast(Sv,1-exp(-P.a*P.bhv*Ih));

if newVE>0 && sum(w)>0
    pos=find(w>0);
    cdf=cumsum(w(pos)); cdf=cdf/cdf(end);
    bins=discretize(rand(newVE,1),[0;cdf]);
    bins(isnan(bins))=numel(pos);
    donors=pos(bins);
    E.Zhv=accumarray(donors,1,[N 1]);
end

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

function y=logitSafe(p)
p=min(max(p,1e-8),1-1e-8);
y=log(p./(1-p));
end

function p=invlogit(x)
p=1./(1+exp(-x));
end

function g=groupIndex(sex,age)
% 1=fy, 2=fa, 3=my, 4=ma
g=4*ones(numel(age),1);
g(sex & age<5)=1;
g(sex & age>=5)=2;
g(~sex & age<5)=3;
g(~sex & age>=5)=4;
end

function dotWhisker(ax,Y,names,ylab)
[R,K]=size(Y);
for k=1:K
    y=Y(:,k);
    rng(88000+k,'twister'); % plotting jitter only
    xj=k+0.11*(rand(R,1)-0.5);

    scatter(ax,xj,y,13,'filled',...
        'MarkerFaceColor',[0.45 0.60 0.72],...
        'MarkerFaceAlpha',0.10,...
        'MarkerEdgeColor','none',...
        'HandleVisibility','off');

    m=mean(y,'omitnan');
    se=std(y,'omitnan')/sqrt(sum(isfinite(y)));
    lo=m-1.96*se; hi=m+1.96*se;

    errorbar(ax,k,m,m-lo,hi-m,'o',...
        'Color',[0.05 0.05 0.05],...
        'MarkerFaceColor',[0.20 0.50 0.75],...
        'MarkerEdgeColor',[0.05 0.05 0.05],...
        'MarkerSize',6,'LineWidth',1.6,'CapSize',8);
end
xticks(ax,1:K);
xticklabels(ax,names);
xtickangle(ax,18);
ylabel(ax,ylab);
end

function formatPanel(ax)
grid(ax,'on');
ax.GridAlpha=0.15;
ax.LineWidth=0.8;
ax.TickDir='out';
box(ax,'off');
end

function q=localPercentile(x,p)
x=sort(x(:));
x=x(isfinite(x));
if isempty(x), q=NaN; return; end
r=1+(numel(x)-1)*p/100;
lo=floor(r); hi=ceil(r);
if lo==hi
    q=x(lo);
else
    q=x(lo)+(r-lo)*(x(hi)-x(lo));
end
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
    L=exp(-lambda); k=0; p=1;
    while p>L, k=k+1; p=p*rand; end
    x=k-1;
else
    x=max(0,round(lambda+sqrt(lambda)*randn));
end
end
