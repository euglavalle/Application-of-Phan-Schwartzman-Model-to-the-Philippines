delete(gcp('nocreate'));
clear all; clc; close all;
%parpool('local', 8);
%% Calibration
% Basic parameters
%scaleMEX = 50.12724, shapeMEX = 3.054898
%par.costy = 0.075; %.1; % Default cost constant for Bakkensen & Barrage marginal damage mps
%par.psi = 1; % Default cost curvature for Bakkensen & Barrage marginal damage mps

%scalePHL = 48.9662, shapePHL = 5.11328
%par.costy = 0.08;%.1; % Default cost constant
%par.psi = 8.6; % Default cost curvature
par.costy = 0.08;%.1; % Default cost constant for Bakkensen & Barrage marginal damage mps
par.psi = 8.5; % Default cost curvature for Bakkensen & Barrage marginal damage mps

%scalePHL = 95.18268, shapePHL = 5.11328
%par.costy = 0.075;%.1; % Default cost constant for Bakkensen & Barrage marginal damage kts
%par.psi = 9.1; % Default cost curvature for Bakkensen & Barrage marginal damage kts

%scalePHL = 17.246, shapePHL = .992675
%par.costy = 0.08; % Default cost constant
%par.psi = 9; % Default cost curvature

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
%mug = .0066;
%sdg = .0213;

% TFP process: based on quarterly AG JPE 2007 Table 4 column 4 of sdg=.0213
% Note: these are quarterly parameters, and will be adjusted below.
mug = .0065; %value for Philippines from EM-DAT estimation
sdg = .0236; %value for Philippines from EM-DAT estimation
%mug = .0037; %value for Mexico from EM-DAT estimation
%sdg = .0224; %value for Mexico from EM-DAT estimation
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

% Clean unused variables from memory
clear Gammas Gammas_beforelast Gammas_last epsilons gs z

% Build grids
G = 100; B=100; K=100; C=1; M=1000; CC = 50; Z =1;

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
%par.marginaldamage = ((0.0895/100)/par.alpha)*(par.period/5);%marginal GDP damage of 1m/s is 0.0895% of GDP over 5 years
par.marginaldamage = (((0.46/70)/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/70mps, where 0.46% is TFP damage over 5 years per Cat5 storm
%par.marginaldamage = (((0.46/137)/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/137kts, where 0.46% is TFP damage over 5 years per Cat5 storm
%support = [0:50]; %orginal support in mps from Phan & Schwartzman for MEX
%support = [0:100]; %support for PHL. Bigger range than MEX to get rid of spikes in PDF of increased cyclone risk scenarioe
% Define support according to simulated annual maximum wind speed
% values
file_path = 'Weibull_estimation_Barrage/PHL.csv';
%file_path = 'Weibull_estimation_Barrage/MEX.csv';
data = readtable(file_path);
wind_speeds = data.WSms;
%wind_speeds = data.WSkts;
support = linspace(min(wind_speeds), max(wind_speeds), 100);
shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage for mps wind speeds
scale = 48.9662; %current value for Philippines based on Bakkensen & Barrage for mps wind speeds
%shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
%scale = 95.18268; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
%shape = 3.054898; %current value for Mexico based on Bakkensen & Barrage for mps wind speeds 
%scale = 50.12724; %current value for Mexico based on Bakkensen & Barrage for mps wind speeds 
%shape = .992675;% shape and scale chosen to match Hsiang Jina's stats
                    % for windspeed of 90 and 99pctile cyclones (19.5m/s and 39.2m/s)
%shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage
%scale = 95.18268; %current value for Philippines based on Bakkensen & Barrage
%shape = 3.012625; %current value for Mexico based on Bakkensen & Barrage
%scale = 97.26602; %current value for Mexico based on Bakkensen & Barrage
%shape = 6.1008945; %future value for Philippines based on Bakkensen & Barrage
%scale = 108.21446; %future value for Philippines based on Bakkensen & Barrage
%shape = 5.7547312; %future value for Mexico based on Bakkensen & Barrage
%scale = 103.06419; %future value for Mexico based on Bakkensen & Barrage
%scale = 17.246; %value for Philippines based on Hsiang & Jina
%scale = 6.57933;% to match current area-normalized ws of 6.6m/s for Mexico;%8.41678;
par.pperyear = 0.58; % prob of a cyclone landfall per year
                     % =5.8/10 (Hsiang Jina prob 90pct cyclone is 5.8%)
par.p = 1-(1-par.pperyear)^par.period;
par.f = 1;%duration of a short period of time in which we assume at most one cyclone can make landfall
      %f=1/365 is a day; f=1 is a year
%dgrid=disweibull(shape,scale,support,par,10);
dgrid=disweibull_nostorm_included(shape,scale,support,par,10);
dgrid_nonlin = dgrid;
dgrid_nonlin.d(1:round(length(dgrid_nonlin.d)/2)-1)=0;
%% %%%%%%%%%%%%%%%%%%%%%%%%% Section 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Baseline model
out.base = opt_w(par,grid,dgrid,0,ones(M,1));

% Comput baseline egordic steady state
%rng(0); %fix a seed

% Ensure reproducibility in parallel execution
%spmd
    %stream = RandStream.create('mrg32k3a', 'NumStreams', numlabs, 'StreamIndices', labindex);
    %RandStream.setGlobalStream(stream);
%end

Tinf=10^5; par.Tinf=Tinf;
mi = 1; 
tic
sim1g=randgt(grid.pdfg,1,Tinf);
sim1z=randgt(grid.pdfz,1,Tinf);
sim1zval = grid.z(sim1z);

sim1x=dgrid.d(randgt(dgrid.pdf,1,Tinf));
foo1=simulate(sim1g,sim1x,sim1zval,mi,par,grid,dgrid,out.base)
toc

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

out.base

% %%%%%%%%%%%%%%%%%%%%%% Compute IRF to cyclone activity shock %%%%%%%%
N = 10^6;%Number of simulations 
T = 100;%40;%Length of each simulation
tshock=round(2/3*T); %Period of cyclone shock

% Generate random shocks % This will take a few minutes
rng(1); %fix a seed
tic
simg=randgt(grid.pdfg,N,T);
simz=randgt(grid.pdfz,N,T);
simzval = grid.z(simz);
simx=dgrid.d(randgt(dgrid.pdf,N,T));
toc
% Define impulse
par.xshock = nanstd(simx(:,end)); %one std shock 
simx(:,tshock)=simx(:,tshock)+par.xshock;

% Set initial net worth
mi = snapmat(out.base.avem,grid.m); % at egordic average
% Sim
simm=nan(N,T);simk=nan(N,T);simby=nan(N,T);
simnx=nan(N,T);sims=nan(N,T);simdef=nan(N,T);simc=nan(N,T);
disp('Simulation is running, this may take 15 minutes...')
tic
parfor i=1:N
    if mod(i,round(N/1000))==0
    disp(i);
    end
    foo=simulate(simg(i,:),simx(i,:),simzval(i,:),mi,par,grid,dgrid,out.base);
    simm(i,:)=foo.m; simk(i,:)=foo.k; simby(i,:)=foo.by;
    simnx(i,:)=foo.nx;sims(i,:)=foo.s;simdef(i,:)=foo.def;
    simc(i,:)=foo.c;
end
toc

% After parallel computations, reset to sequential RNG (if needed)
%delete(gcp('nocreate')); % Close parallel pool if not needed anymore
%rng(1); % Reset for sequential calculations

% Generating Figure 3: IRF to cyclone shock
tmin=tshock-1; tmax=tshock+10;
tt = par.period*(tmin-tshock:tmax-tshock);

out.irfx = mean(simx(:,tmin:tmax),1)/par.marginaldamage;
out.irfk = mean(simk(:,tmin:tmax),1); out.irfk = 100*(out.irfk./out.base.avek-1);
out.irfy = mean(simk(:,tmin:tmax).^par.alpha,1); out.irfy = 100*(out.irfy./out.base.avey-1);
out.irfby = 100*par.period*mean(simby(:,tmin:tmax),1);
out.irfs = 100/par.period*mean(sims(:,tmin:tmax),1);
out.irfdef = 100/par.period*mean(simdef(:,tmin:tmax),1);
out.irfnx = -100*par.period*mean(simnx(:,tmin:tmax),1);
out.irfm = mean(simm(:,tmin:tmax),1); out.irfm = 100*(out.irfm./out.base.avem-1);
out.irfc = mean(simc(:,tmin:tmax),1); out.irfc = 100*(out.irfc./out.base.avec-1);

%disp('Printing out.irfby and out.irfs')
%out.irfby
%out.irfs

% 
figure
row = 2; col = 3; 
subplot(row,col,1); 
plot(tt,out.irfx);title('Cyclone activity'); 
axis([-inf,inf,0,inf]); xlabel('years'); %ylabel('knots');
ylabel('m/s');

subplot(row,col,2); 
plot(tt,out.irfk, tt,out.irfy,'--'); title('Capital & Output')
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% deviation from trend');
text(2,0.95*min(out.irfk),'k');text(0,0.99*min(out.irfy),'y');

subplot(row,col,3);
plot(tt,out.irfm, tt,out.irfc,'--'); title('Net worth & Consumption');
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% deviation from trend');
text(2,0.95*min(out.irfc),'c');text(0,0.6*min(out.irfm),'m');

subplot(row,col,4); 
plot(tt,out.irfby); title('Debt/Output'); 
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% (annualized)');

subplot(row,col,5); 
plot(tt,out.irfs, tt,out.irfdef,'--'); title('Spread & Default');
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% (annualized)');
text(0,0.99*max(out.irfs),'s');text(5,.99*max(out.irfdef),'def');

subplot(row,col,6);
plot(tt,out.irfnx); title('Capital inflows/Output'); 
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% (annualized)');

saveas(gcf,'output_phl\scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm\irfstd.png');
%% %%%%%%%%%%%%%%%%%%%%%% Section 2 Climate Change %%%%%%%%%%%%%%%%%%%%%%%%
% Generate Fig 4

% Climate change leads to 10\% increase in cyclone activity
%dgridcc=disweibull(shape,scale*1.10,support,par,10); 
%dgridcc=disweibull(shape,scale*1.191,support,par,10); %increase in western Pacific basin according to Emanuel et al., 2008
%dgridcc=disweibull_nostorm_included(shape,scale*1.10,support,par,10);
dgridcc=disweibull_nostorm_included(shape,scale*1.191,support,par,10); %increase in western Pacific basin according to Emanuel et al., 2008
%dgridcc=disweibull_nostorm_included(6.1008945,108.21446,support,par,10); %future shape & scale parameters for PHL
%dgridcc=disweibull_nostorm_included(6.5302267,95.69236,support,par,10); %future shape & scale parameters for PHL ECHAM
%dgridcc=disweibull_nostorm_included(5.9966989,95.371552,support,par,10); %future shape & scale parameters for PHL CNRM
%dgridcc=disweibull_nostorm_included(5.0599093,108.89911,support,par,10); %future shape & scale parameters for PHL MIROC

dgridcc_notail = dgridcc;
dgridcc_notail.d(end-1:end) = dgrid.d(end-1:end); %replace last two damage states with values from dgrid.d
dgridcc_notail.pdf = dgridcc.pdf; %pdf unchanged, i.e. probability assigned to each state in dgridcc_notail 
% is the same as in dgridcc. -> preserves the overall distribution shape and ensures consistency in probability weights.
dgridcc_notail.d = dgridcc_notail.d.*sum(dgridcc.d.*dgridcc.pdf)./sum(dgridcc_notail.d.*dgridcc_notail.pdf);
%scaling factor adjusts dgridcc_notail.d such that the overall mean damage remains the same before and after modifying the tail.

% Plot Fig 4a for PDF and CDF of cyclone activitity after climate change
figure
row=1;col=2;
subplot(row,col,1)
plot(dgrid.dwspeed,dgrid.pdfwspeed,'b--',...
    dgridcc.dwspeed,dgridcc.pdfwspeed,'r-','LineWidth',2);
%xlabel('Cyclone activity (5-year cumulative, knots)');
xlabel('Cyclone activity (5-year cumulative, m/s)'); 
title('PDF');

subplot(row,col,2)
plot(dgrid.dwspeed,dgrid.cdfwspeed,'b--', ...
    dgridcc.dwspeed,dgridcc.cdfwspeed,'r-','LineWidth',2);
xlabel('Cyclone activity'); 
title('CDF');
legend('Baseline','Increased cyclone risk','Location','best');

saveas(gcf,'output_phl\scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm\pdf.png');

% Climate change welfare loss
out.cc = opt_w(par,grid,dgridcc,0,out.base.v); %take out.base.v as initial guess for computational efficiency
loss = 100*(1-out.cc.v./out.base.v);

% Plot Figure 4b: histograms and welfare functions
rng(0); %fix a seed
simx1cc=dgridcc.d(randgt(dgridcc.pdf,1,Tinf));
foocc1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.cc)

data1 = foo1.m;
data2 = foocc1.m;
% Smooth the data using ksdensity
[f1, x1] = ksdensity(data1);
[f2, x2] = ksdensity(data2);

figure;
hold on;
% Plot welfare functions
yyaxis left;
plot(grid.m,out.base.v,'b--',...
    grid.m,out.cc.v,'r-','LineWidth',2);
ylabel('Welfare function v');
xmin = out.base.avem - 3*out.base.stdm;
xmax = out.base.avem + 3*out.base.stdm;
ymin = out.base.v(snapmat(xmin,grid.m));
ymax = out.base.v(snapmat(xmax,grid.m));
ylim([ymin,ymax]);
xlim([xmin,xmax]);

% Create the histograms with more bins and normalization
transparency = 0.1;
yyaxis right;
fill([x1, fliplr(x1)], [f1, zeros(size(f1))], 'b', 'EdgeColor','b',...
    'FaceAlpha', transparency,'LineStyle', '--');
fill([x2, fliplr(x2)], [f2, zeros(size(f2))], 'r', 'EdgeColor','r',...
    'FaceAlpha', transparency,'LineStyle', '-');
ylabel('Density');
ylabel('Frequency'); 
%ylim([0,5]);yticks([0,1,2,3]);
hold off

% Set the axis labels and legend
xlabel('Net worth m');
ylabel('Frequency');
legendEntries={...%sprintf('Welfare loss 1-v_+(m)/v(m)'),...%
     sprintf('v(m)'),sprintf('v_+(m)'),...
    sprintf('PDF(m)'),...
    sprintf('PDF_+(m)')};
legend(legendEntries,... 
    'Location', 'northoutside', 'Orientation', 'horizontal',...
    'NumColumns', 2, 'Interpreter', 'tex');

%Add titles to legend
[leg,att]=legend('show');
title(leg,sprintf('Baseline \t            Increased cyclone risk'));
leg.Title.Visible='on';

saveas(gcf,'output_phl\scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm\v_hist.png');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Section 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Financial adaptation

% Insurance in climate changed world
out.ins = opt_w(par,grid,dgridcc,1,out.cc.v)
gain.ins = 100*(out.ins.v./out.cc.v-1);

% CAT bond + disaster insurance in climate changed world
% Note: Solving this optimization takes >1 hour on 24 workers
% on Intel(R) Xeon(R) Silver 4214 CPU @ 2.20GHz   2.19 GHz  (2 processors)
dgridcc.catthreshold = dgridcc.d(find(dgridcc.cdfcond > .90,1))%cutoff is 90% pctile of nonzero cyclone damage
                %0;%any cyclone would trigger cat clause
gridcat=grid;
maxcat=0.25; C=CC; 
gridcat.cat = [0:maxcat/(C-1):maxcat]';

out.catins = opt_w(par,gridcat,dgridcc,1,out.cc.v)
gain.catins = 100*(out.catins.v./out.cc.v-1);

% CAT bond only
% Note: Solving this optimization takes >1 hour on 24 workers
out.cat = opt_w(par,gridcat,dgridcc,0,out.catins.v)
gain.cat = 100*(out.cat.v./out.cc.v-1);

% Plot welfare gain from financial adaptation
loss = 100*(1-out.cc.v./out.base.v);

figure
plot(grid.m, smoothdata(loss),'r',...
    grid.m, smoothdata(gain.catins),'--g',...
    grid.m, smoothdata(gain.ins),'-.b',...
    grid.m, smoothdata(gain.cat),':k');
xlabel('m'); ylabel('% permanent consumption');
title({'Welfare gain from financial adaptation'}); 
legend('Loss from increased cyclone risk',...
        'Gain from CAT bonds + insurance',...
      'Gain from disaster insurance',...
      'Gain from CAT',  'location','best');
xlim([0.7,1.6]);

saveas(gcf,'output_phl\scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm\v_catins.png');

%  ERGODIC WELFARE COMPARISON
rng(0); %fix a seed
simx1cc=dgridcc.d(randgt(dgridcc.pdf,1,Tinf));

mi = 1; 
% Simulate 
tic
foocc1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.cc)
foocc_notail1 = simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc_notail,out.cc)
fooins1=simulate(sim1g,simx1cc,sim1zval,mi,par,grid,dgridcc,out.ins)
foocat1=simulate(sim1g,simx1cc,sim1zval,mi,par,gridcat,dgridcc,out.cat)
foocatins1=simulate(sim1g,simx1cc,sim1zval,mi,par,gridcat,dgridcc,out.catins)

toc
% Welfare comparisons

% averages
out.cc.avem=nanmean(foocc1.m(100:end));
out.cc.aveby=nanmean(foocc1.by(100:end));
out.cc.avek=nanmean(foocc1.k(100:end));
out.cc.avey=nanmean(foocc1.k(100:end).^par.alpha);
out.cc.aves = nanmean(foocc1.s(100:end));
out.cc.avedef = nanmean(foocc1.def(100:end));
out.cc.avev = nanmean(foocc1.v(100:end));

out.cc_notail.avem=nanmean(foocc_notail1.m(100:end));
out.cc_notail.aveby=nanmean(foocc_notail1.by(100:end));
out.cc_notail.avek=nanmean(foocc_notail1.k(100:end));
out.cc_notail.avey=nanmean(foocc_notail1.k(100:end).^par.alpha);
out.cc_notail.aves = nanmean(foocc_notail1.s(100:end));
out.cc_notail.avedef = nanmean(foocc_notail1.def(100:end));
out.cc_notail.avev = nanmean(foocc_notail1.v(100:end));

out.ins.avem=nanmean(fooins1.m(100:end));
out.ins.aveby=nanmean(fooins1.by(100:end));
out.ins.avek=nanmean(fooins1.k(100:end));
out.ins.avey=nanmean(fooins1.k(100:end).^par.alpha);
out.ins.aves = nanmean(fooins1.s(100:end));
out.ins.avedef = nanmean(fooins1.def(100:end));
out.ins.avev = nanmean(fooins1.v(100:end));

out.cat.avem=nanmean(foocat1.m(100:end));
out.cat.aveby=nanmean(foocat1.by(100:end));
out.cat.avek=nanmean(foocat1.k(100:end));
out.cat.avey=nanmean(foocat1.k(100:end).^par.alpha);
out.cat.aves = nanmean(foocat1.s(100:end));
out.cat.avedef = nanmean(foocat1.def(100:end));
out.cat.avev = nanmean(foocat1.v(100:end));

out.catins.avem=nanmean(foocatins1.m(100:end));
out.catins.aveby=nanmean(foocatins1.by(100:end));
out.catins.avek=nanmean(foocatins1.k(100:end));
out.catins.avey=nanmean(foocatins1.k(100:end).^par.alpha);
out.catins.aves = nanmean(foocatins1.s(100:end));
out.catins.avedef = nanmean(foocatins1.def(100:end));
out.catins.avev = nanmean(foocatins1.v(100:end));

% std
out.cc.stdm=nanstd(foocc1.m(100:end));
out.cc.stdby=nanstd(foocc1.by(100:end));
out.cc.stdk=nanstd(foocc1.k(100:end));
out.cc.stdy=nanstd(foocc1.k(100:end).^par.alpha);
out.cc.stds = nanstd(foocc1.s(100:end));
out.cc.stddef = nanstd(foocc1.def(100:end));
out.cc.stdv = nanstd(foocc1.v(100:end));

out.ins.stdm=nanstd(fooins1.m(100:end));
out.ins.stdby=nanstd(fooins1.by(100:end));
out.ins.stdk=nanstd(fooins1.k(100:end));
out.ins.stdy=nanstd(fooins1.k(100:end).^par.alpha);
out.ins.stds = nanstd(fooins1.s(100:end));
out.ins.stddef = nanstd(fooins1.def(100:end));
out.ins.stdv = nanstd(fooins1.v(100:end));

out.cat.stdm=nanstd(foocat1.m(100:end));
out.cat.stdby=nanstd(foocat1.by(100:end));
out.cat.stdk=nanstd(foocat1.k(100:end));
out.cat.stdy=nanstd(foocat1.k(100:end).^par.alpha);
out.cat.stds = nanstd(foocat1.s(100:end));
out.cat.stddef = nanstd(foocat1.def(100:end));
out.cat.stdv = nanstd(foocat1.v(100:end));

out.catins.stdm=nanstd(foocatins1.m(100:end));
out.catins.stdby=nanstd(foocatins1.by(100:end));
out.catins.stdk=nanstd(foocatins1.k(100:end));
out.catins.stdy=nanstd(foocatins1.k(100:end).^par.alpha);
out.catins.stds = nanstd(foocatins1.s(100:end));
out.catins.stddef = nanstd(foocatins1.def(100:end));
out.catins.stdv = nanstd(foocatins1.v(100:end));

welfaretab = 100*[
out.cc.avem./out.base.avem-1,  out.cc.avev./out.base.avev-1;
%out.notail.avem./out.cc.avem-1, out.notail.avev./out.cc.avev-1;
out.ins.avem./out.cc.avem-1,  out.ins.avev./out.cc.avev-1;
out.cat.avem./out.cc.avem-1,  out.cat.avev./out.cc.avev-1;
out.catins.avem./out.cc.avem-1,  out.catins.avev./out.cc.avev-1;
]

folder_path = 'output_phl/scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm';
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

file_path = fullfile(folder_path, 'welfaretab.csv');
writematrix(welfaretab,file_path);
%% %%%%%%%%%%%%%%%%%%%%%%%%% Appendix Figure A1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cyclone only iid shock to tfp
%Redo grids - take out d shock, put in z shock
% Cyclone strength follows Weibull distribution
% Disaster parameters from Hsiang Jina

par_z = par; %Hurricane distribution does not change, only marginal damage to z vs k
%par_z.marginaldamage = (0.0895/100)*(par.period/5);%marginal GDP damage of 1m/s is 0.0895% of GDP over 5 years. 
par.marginaldamage = (((0.46/70)/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/70mps, where 0.46% is TFP damage over 5 years per Cat5 storm
%par_z.marginaldamage = (((0.46/137)/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/137kts, where 0.46% is TFP damage over 5 years per Cat5 storm
% Damage to k is scaled by 1/alpha to match gdp. Damage to z is not scaled

%gridz_z = disweibull(shape,scale,support,par_z,10);
gridz_z=disweibull_nostorm_included(shape,scale,support,par_z,10);
grid_z = grid;
grid_z.z = -gridz_z.d';
grid_z.pdfz = gridz_z.pdf;
grid_z.pdfz = grid_z.pdfz(1,:);
cdfz_z = cumsum(grid_z.pdfz);
%Now there is no d shock, only z shock
dgrid_z.d=0;
dgrid_z.pdf = 1;

%%% Baseline model
out.base_z = opt_w(par_z,grid_z,dgrid_z,0,ones(M,1));

% Comput baseline egordic steady state
rng(0); %fix a seed
Tinf=10^5; par.Tinf=Tinf;
mi = 1; 
tic
sim1g_z=randgt(grid_z.pdfg,1,Tinf);
sim1z_z=randgt(grid_z.pdfz,1,Tinf);
sim1zval_z = grid_z.z(sim1z_z)';
sim1x_z=dgrid.d(randgt(dgrid_z.pdf,1,Tinf));
foo1=simulate(sim1g_z,sim1x_z,sim1zval_z,mi,par_z,grid_z,dgrid_z,out.base_z)
toc

% averages
out.base_z.avem=nanmean(foo1.m(100:end));
out.base_z.aveby=nanmean(foo1.by(100:end));
out.base_z.avek=nanmean(foo1.k(100:end));
out.base_z.avey=nanmean(exp(sim1zval_z(100:end)).*foo1.k(100:end).^par.alpha);
out.base_z.aves = nanmean(foo1.s(100:end));
out.base_z.avedef = nanmean(foo1.def(100:end));
out.base_z.avev = nanmean(foo1.v(100:end));
out.base_z.avec = nanmean(foo1.c(100:end));

out.base_z.stdm=nanstd(foo1.m(100:end));
out.base_z.stdby=nanstd(foo1.by(100:end));
out.base_z.stdk=nanstd(foo1.k(100:end));
out.base_z.stdy=nanstd(exp(sim1zval_z(100:end)).*foo1.k(100:end).^par.alpha);
out.base_z.stds = nanstd(foo1.s(100:end));
out.base_z.stddef = nanstd(foo1.def(100:end));
out.base_z.stdv = nanstd(foo1.v(100:end));
out.base_z.stdc = nanstd(foo1.c(100:end));

out.base_z

%%% IRF to cyclone activity shock 
N = 10^6;
T = 100;
tshock=round(2/3*T); %Period of cyclone shock

% Generate random shocks % This will take 1.5 minutes
rng(1); %fix a seed
tic
simg_z=randgt(grid_z.pdfg,N,T);
simz_z=randgt(grid_z.pdfz,N,T);
simx_z=dgrid_z.d(randgt(dgrid_z.pdf,N,T));
toc
% Define impulse
simzval_z = grid_z.z(simz_z);
par_z.xshock = nanstd(simzval_z(:,end)); %one std shock
simzval_z(:,tshock)=simzval_z(:,tshock)-par_z.xshock;

% Set initial net worth
mi = snapmat(out.base_z.avem,grid_z.m); % at egordic average
% Sim
simm_z=nan(N,T);simk_z=nan(N,T);simby_z=nan(N,T);
simnx_z=nan(N,T);sims_z=nan(N,T);simdef_z=nan(N,T);simc_z=nan(N,T);
tic
parfor i=1:N
    if mod(i,round(N/1000))==0
    disp(i);
    end
    foo=simulate(simg_z(i,:),simx_z(i,:),simzval_z(i,:),mi,par_z,grid_z,dgrid_z,out.base_z);
    simm_z(i,:)=foo.m; simk_z(i,:)=foo.k; simby_z(i,:)=foo.by;
    simnx_z(i,:)=foo.nx;sims_z(i,:)=foo.s;simdef_z(i,:)=foo.def;
    simc_z(i,:)=foo.c;
end
toc

tmin=tshock-1; tmax=tshock+10;
tt = par.period*(tmin-tshock:tmax-tshock);

out.irfx_z = mean(simzval_z(:,tmin:tmax),1);
out.irfz_z = mean(simzval_z(:,tmin:tmax),1);
out.irfk_z = mean(simk_z(:,tmin:tmax),1); out.irfk_z = 100*(out.irfk_z./out.base_z.avek-1);
out.irfy_z = mean(exp(simzval_z(:,tmin:tmax)).*simk_z(:,tmin:tmax).^par.alpha,1); out.irfy_z = 100*(out.irfy_z./out.base_z.avey-1);
out.irfby_z = 100*par.period*mean(simby_z(:,tmin:tmax),1);
out.irfs_z = 100/par.period*mean(sims_z(:,tmin:tmax),1);
out.irfdef_z = 100/par.period*mean(simdef_z(:,tmin:tmax),1);
out.irfnx_z = -100*par.period*mean(simnx_z(:,tmin:tmax),1);
out.irfm_z = mean(simm_z(:,tmin:tmax),1); out.irfm_z = 100*(out.irfm_z./out.base_z.avem-1);
out.irfc_z = mean(simc_z(:,tmin:tmax),1); out.irfc_z = 100*(out.irfc_z./out.base_z.avec-1);

% Generate Fig A1
figure
row = 1; col = 4; 
subplot(row,col,1); 
plot(tt,out.irfy, tt,out.irfy_z,'--'); title('Output')
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% deviation from trend');

subplot(row,col,2); 
plot(tt,out.irfk, tt,out.irfk_z,'--'); title('Capital')
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% deviation from trend');

subplot(row,col,3);
plot(tt,out.irfm, tt,out.irfm_z,'--'); title('Net worth');
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% deviation from trend');

subplot(row,col,4); 
plot(tt,out.irfs, tt,out.irfs_z,'--'); title('Spread');
axis([-inf,inf,-inf,inf]); xlabel('years'); ylabel('% (annualized)');

saveas(gcf,'output_phl\scalecur48_9662_shape5_11328_mugfromEMDAT_newmu_1Cat5storm\irfstd_A1.png');