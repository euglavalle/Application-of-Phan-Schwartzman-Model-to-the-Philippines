function out = runme_in_grid_search(costy, psi)
    %% Calibration
    % Basic parameters
    par.costy = costy;%.1; % Default cost constant
    par.psi = psi; % Default cost curvature
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
    mug = .0065;
    sdg = .0236;
    %mug = .0037;
    %sdg = .0224;
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
    %par.marginaldamage = ((0.0065/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/70mps, where 0.46% is TFP damage over 5 years per Cat5 storm
    par.marginaldamage = (((0.46/137)/100)/par.alpha)*(par.period/5);%use Bakkensen & Barrage estimate of 0.46%/137kts, where 0.46% is TFP damage over 5 years per Cat5 storm
    %support = [0:50];
    % Define support according to simulated annual maximum wind speed
    % values
    file_path = 'Weibull_estimation_Barrage/PHL.csv';
    %file_path = 'Weibull_estimation_Barrage/MEX.csv';
    data = readtable(file_path);
    %wind_speeds = data.WSms;
    wind_speeds = data.WSkts;
    support = linspace(min(wind_speeds), max(wind_speeds), 100);
    %
    %shape = .992675;% shape and scale chosen to match Hsiang Jina's stats
                    % for windspeed of 90 and 99pctile cyclones (19.5m/s and 39.2m/s)
    %shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage for mps wind speeds
    %scale = 48.9662; %current value for Philippines based on Bakkensen &
    %Barrage for mps wind speeds
    shape = 5.11328; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
    scale = 95.18268; %current value for Philippines based on Bakkensen & Barrage for kts wind speeds
    %scale = 95.18268; %current value for Philippines based on Bakkensen & Barrage
    %shape = 3.012625; %current value for Mexico based on Bakkensen & Barrage
    %scale = 97.26602; %current value for Mexico based on Bakkensen & Barrage
    %shape = 3.054898; %current value for Mexico based on Bakkensen & Barrage
    %for mps wind speeds 
    %scale = 50.12724; %current value for Mexico based on Bakkensen & Barrage
    %for mps wind speeds 
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
    spmd
        %stream = RandStream.create('mrg32k3a', 'NumStreams', numlabs, 'StreamIndices', labindex);
        stream = RandStream.create('mrg32k3a', 'NumStreams', spmdSize, 'StreamIndices', spmdIndex);
        RandStream.setGlobalStream(stream);
    end
    
    Tinf=10^5; par.Tinf=Tinf;
    mi = 1; 
    %tic
    sim1g=randgt(grid.pdfg,1,Tinf);
    sim1z=randgt(grid.pdfz,1,Tinf);
    sim1zval = grid.z(sim1z);
    
    sim1x=dgrid.d(randgt(dgrid.pdf,1,Tinf));
    foo1=simulate(sim1g,sim1x,sim1zval,mi,par,grid,dgrid,out.base);
    %toc
    
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
    
    out.base;

    %disp('Fields in out.base:');
    %disp(fieldnames(out.base));