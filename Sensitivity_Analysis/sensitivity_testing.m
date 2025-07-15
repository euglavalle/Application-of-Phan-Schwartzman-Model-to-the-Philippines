delete(gcp('nocreate'));
clear all; clc; close all;

%scalePHL = 17.246, shapePHL = .992675
%par.costy = 0.08; % Default cost constant
%par.psi = 9; % Default cost curvature

%scalePHL = 95.18268, shapePHL = 5.11328
par.costy = 0.08;%.1; % Default cost constant for Bakkensen & Barrage marginal damage kts
par.psi = 8.4; % Default cost curvature for Bakkensen & Barrage marginal damage kts

par.ez=1; %choose EZ specification
par.iies = 0.5; 
par.rra = 4;
par.s = 2;% risk aversion for CRRA specification
par.alpha = 2/3;
par.period = 5; %number of years in a period
par.delta = 1-.9^par.period; %10% depreciation per year
par.r = 1.01^par.period-1; %1% per year
par.beta = .96^par.period; %.96 per year

% TFP process: based on quarterly AG JPE 2007 Table 4 column 4 of sdg=.0213
% Note: these are quarterly parameters, and will be adjusted below.
mug = .0065; %value for Philippines from EM-DAT estimation
sdg = .0236; %value for Philippines from EM-DAT estimation
sdz = .0053;
rhoz = .95;

% Simulate cumulative shocks Gamma 
rng(0)
gs = mug+sdg*randn(1e6,100);
Gammas = ones(size(gs));
for t = 2:100
    Gammas(:,t) = Gammas(:,t-1).*exp(gs(:,t));
end

% Time-aggregate Gammas from a quarter to a period of 5 years.
% To do this, think about output in a given 5-year period as being the sum
% of output in all the quarters within that period. For our time aggregation we
% calculate what output would be in each of those periods for fixed k and z.
% This means that we can define time aggregated Gamma_bar to b proportional to 
% \sum_{t in time aggregation}Gamma_{t}, and then back out the process for
% Gamma_bar from the Montecarlo. We do the the same for z shocks later.
Gammas_last = sum(Gammas(:,end-par.period*4+1:end),2);
Gammas_beforelast = sum(Gammas(:,end-2*par.period*4+1:end-par.period*4),2);
par.mug = mean(Gammas_last./Gammas_beforelast);
par.sdg = std(Gammas_last./Gammas_beforelast);
par.rhog = 0;

% Time-aggregate transitory shocks Zs
epsilons = sdz*randn(1e6,100);
z = zeros(size(epsilons));
for t=2:size(epsilons,2)
z(:,t)=rhoz*z(:,t-1)+epsilons(:,t);
end
Z = log(sum(exp(z(:,t-19:end)),2));
par.sdz = std(Z);
par.tolv=1e-6;
%par.tolv=1e-8;

% Clean unused variables from memory
clear Gammas Gammas_beforelast Gammas_last epsilons gs z

% Build grids
G = 100; B=100; K=100; C=1; M=3000; CC = 50; Z =1;

kmin=0.01;
kmax = 2;
[grid.g,grid.pdfg]=dis(par.mug, par.sdg, par.rhog, G); 
[grid.z,grid.pdfz]=disz(0, par.sdz, 0, Z); 
grid.pdfg=grid.pdfg(1,:);%iid shocks;
grid.pdfz = grid.pdfz(1,:);
cdfg=cumsum(grid.pdfg);
cdfz = cumsum(grid.pdfz);

bymin= par.costy*min(grid.g)/exp(max(grid.z)); 
bymax = par.costy*max(grid.g)/exp(min(grid.z));
grid.by=[0,[bymin:(bymax-bymin)/(B-2):bymax]]';
grid.k=[kmin:(kmax-kmin)/(K-1):kmax]';
grid.b = grid.by.*median(grid.k).^par.alpha;

if C==1
    grid.cat=[0];
else
    grid.cat = [0:1/(C-1):1]';
end

mmin = kmin;
mmax = (1-par.costy)*(kmax/min(grid.g))^par.alpha+(1-par.delta)*kmax/min(grid.g);
grid.m = [mmin:(mmax-mmin)/(M-1):mmax]';

% Create uniformly sized grids for choices. Ordering: b,cat,k
grid.kn=kron(grid.k,kron(ones(C,1),ones(B,1))); 
grid.byn=kron(ones(K,1),kron(ones(C,1),grid.by)); 
grid.bn = grid.byn.*median(grid.kn).^par.alpha;
grid.catn = kron(ones(K,1),kron(grid.cat,ones(B,1)));

% Cyclone strength follows Weibull distribution
% Disaster parameters from Hsiang Jina
%support = [0:100]; %support for PHL. Bigger range than MEX to get rid of spikes in PDF of increased cyclone risk scenario
%shape = .992675;% shape and scale chosen to match Hsiang Jina's stats
                    % for windspeed of 90 and 99pctile cyclones (19.5m/s and 39.2m/s)
%scale = 17.246; %value for Philippines based on Hsiang & Jina
file_path = 'Weibull_estimation_Barrage/PHL.csv';
data = readtable(file_path);
wind_speeds = data.WSkts;
support = linspace(min(wind_speeds), max(wind_speeds), 100);
shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
scale = 95.18268; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
par.pperyear = 0.58; % prob of a cyclone landfall per year
                     % =5.8/10 (Hsiang Jina prob 90pct cyclone is 5.8%)

par.p = 1-(1-par.pperyear)^par.period;
par.f = 1;%duration of a short period of time in which we assume at most one cyclone can make landfall
      %f=1/365 is a day; f=1 is a year
%par.marginaldamage = ((0.0895/100)/par.alpha)*(par.period/5);%marginal GDP damage of 1m/s is 0.0895% of GDP over 5 years
par.marginaldamage = 0.5*(((0.46/137)/100))*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/137kts, where 0.46% is TFP damage over 5 years per Cat5 storm
%dgrid=disweibull(shape,scale,support,par,.9998,10);
dgrid=disweibull_nostorm_included(shape,scale,support,par,.999,10);


disp('Damage values in dgrid.d:');
disp(dgrid.d);

disp('Prob of values in dgrid.pdf:');
disp(dgrid.pdf);

disp('Min and max of dgrid.d:');
disp([min(dgrid.d), max(dgrid.d)]);
%{
% Check nonzero entries
disp('Number of nonzero damage states:');
disp(sum(dgrid.d > 0));

figure;
bar(dgrid.d, dgrid.pdf);
xlabel('Damage states');
ylabel('Probability');
title('Distribution of Cyclone Damage');
%}

out.base = opt_w(par,grid,dgrid,0,ones(M,1));

Tinf=10^5; par.Tinf=Tinf;
mi = 1; 
tic
sim1g=randgt(grid.pdfg,1,Tinf);
sim1z=randgt(grid.pdfz,1,Tinf);
sim1zval = grid.z(sim1z);

sim1x=dgrid.d(randgt(dgrid.pdf,1,Tinf));
foo1=simulate(sim1g,sim1x,sim1zval,mi,par,grid,dgrid,out.base);
toc

disp('Simulated damage summary:');
fprintf('Min: %.4e, Max: %.4e, Mean: %.4e, Median: %.4e, Mode: %.4e\n', min(sim1x), max(sim1x), mean(sim1x), median(sim1x), mode(sim1x));

% Averages
out.base.avem=nanmean(foo1.m(100:end));
out.base.aveby=nanmean(foo1.by(100:end));
out.base.avek=nanmean(foo1.k(100:end));
out.base.avey=nanmean(foo1.k(100:end).^par.alpha);
out.base.aves = nanmean(foo1.s(100:end));
out.base.avedef = nanmean(foo1.def(100:end));
out.base.avev = nanmean(foo1.v(100:end));
out.base.avec = nanmean(foo1.c(100:end));

out.base.stdm=nanstd(foo1.m(100:end));
out.base.stdby=nanstd(foo1.by(100:end));
out.base.stdk=nanstd(foo1.k(100:end));
out.base.stdy=nanstd(foo1.k(100:end).^par.alpha);
out.base.stds = nanstd(foo1.s(100:end));
out.base.stddef = nanstd(foo1.def(100:end));
out.base.stdv = nanstd(foo1.v(100:end));
out.base.stdc = nanstd(foo1.c(100:end));

rng(1);
%dgridcc=disweibull(shape,scale*1.191,support,par,.9998,10); %increase in western Pacific basin according to Emanuel et al., 2008
dgridcc=disweibull_nostorm_included(shape,scale*1.191,support,par,.999,10); %increase in western Pacific basin according to Emanuel et al., 2008
%dgridcc=disweibull_nostorm_included(5.9966989,95.371552,support,par,.999,10); %future shape & scale parameters for PHL CNRM
%dgridcc=disweibull_nostorm_included(6.5302267,95.69236,support,par,.999,10); %future shape & scale parameters for PHL ECHAM
%dgridcc=disweibull_nostorm_included(6.1008945,108.21446,support,par,.999,10);%future shape & scale parameters for PHL GFDL
%dgridcc=disweibull_nostorm_included(5.0599093,108.89911,support,par,.999,10); %future shape & scale parameters for PHL MIROC

dgridcc_notail = dgridcc;

dgridcc_notail.d(end-1:end) = dgrid.d(end-1:end); %replace last two damage states with values from dgrid.d
dgridcc_notail.pdf = dgridcc.pdf; %pdf unchanged, i.e. probability assigned to each state in dgridcc_notail 
% is the same as in dgridcc. -> preserves the overall distribution shape and ensures consistency in probability weights.
dgridcc_notail.d = dgridcc_notail.d.*sum(dgridcc.d.*dgridcc.pdf)./sum(dgridcc_notail.d.*dgridcc_notail.pdf);
%scaling factor adjusts dgridcc_notail.d such that the overall mean damage remains the same before and after modifying the tail.

% Climate change welfare loss
out.cc = opt_w(par,grid,dgridcc,0,out.base.v); %take out.base.v as initial guess for computational efficiency
loss = 100*(1-out.cc.v./out.base.v);

% Plot Figure 4b: histograms and welfare functions
rng(0); %fix a seed
simx1cc=dgridcc.d(randgt(dgridcc.pdf,1,Tinf));

disp('Simulated climate change damage summary:');
fprintf('Min: %.4e, Max: %.4e, Mean: %.4e, Median: %.4e, Mode: %.4e\n', min(simx1cc), max(simx1cc), mean(simx1cc), median(simx1cc), mode(simx1cc));

foocc1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.cc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Insurance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Financial adaptation

% Insurance in climate changed world
out.ins = opt_w(par,grid,dgridcc,1,out.cc.v);
gain.ins = 100*(out.ins.v./out.cc.v-1);

% Plot welfare gain from financial adaptation
loss = 100*(1-out.cc.v./out.base.v);

%  ERGODIC WELFARE COMPARISON
rng(0); %fix a seed
simx1cc=dgridcc.d(randgt(dgridcc.pdf,1,Tinf));

mi = 1; 
% Simulate 
tic
foocc1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.cc);
foocc_notail1 = simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc_notail,out.cc);
fooins1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.ins);
toc

 % Welfare comparisons
    
% averages
out.cc.avev = nanmean(foocc1.v(100:end));

out.ins.avev = nanmean(fooins1.v(100:end));

fprintf('Welfare gain with insurance: %.6f\n', 100*(out.ins.avev./out.cc.avev-1));