function out=simulate(gt,xt,zt,minit,par,grid,dgrid,inputpolicy)
T=size(gt,2);
NN = size(gt,1);
[B,~] = size(grid.by);
[K,~] = size(grid.k);
[C,~] = size(grid.cat);
cat_trigger_count = zeros(NN,1);
cat_trigger_total = zeros(NN,1);

kt=nan(NN,T); %kt(1) = 1;
byt=nan(NN,T); %byt(1) = 1;
def=zeros(NN,T);
st=nan(NN,T);
qt=nan(NN,T);
kn = nan(NN,T); %kn(T) = 1;
byn = nan(NN,T); %byn(T) = 1;
bt = nan(NN,T); %bt(1)=0;
bn = nan(NN,T);
cat = nan(NN,T);
mt = nan(NN,T);
nx = nan(NN,T);
ct = nan(NN,T);
vt = nan(NN,T);
mt(:,1)=grid.m(minit);
for t=1:T-1
    mi=snapmat(mt(:,t),grid.m);%snap mt to mgrid
    vt(:,t)=inputpolicy.v(mi);
    byn(:,t)=inputpolicy.policybyn(mi);
    bn(:,t)=inputpolicy.policybn(mi);
    cat(:,t)=inputpolicy.policycat(mi);
    kn(:,t)=inputpolicy.policykn(mi);
    bn(:,t)=byn(:,t).*kn(:,t).^par.alpha;
    
    spread = reshape(inputpolicy.sp,B,C,K);

    if NN>1
        parfor n = 1:NN
        st(n,t)=spread(snapmat(byn(n,t),grid.by),...%grids.by==byn(t),...
                       snapmat(cat(n,t),grid.cat),...%grids.cat==cat(t),...
                       snapmat(kn(n,t),grid.k));%grids.k==kn(t));
        end
    else
        st(1,t)=spread(snapmat(byn(1,t),grid.by),...%grids.by==byn(t),...
                       snapmat(cat(1,t),grid.cat),...%grids.cat==cat(t),...
                       snapmat(kn(1,t),grid.k));%grids.k==kn(t));
    end
    qt(:,t)=(1-st(:,t))/(1+par.r);
    
    kt(:,t+1)=exp( -xt(:,t+1) ).*kn(:,t)./grid.g(gt(:,t+1));%k':=exp(-d'*dk)*kn/G'
    
    if C>1
        bt(:,t+1)=( 1-cat(:,t)*(xt(:,t+1)>dgrid.catthreshold) ).*bn(:,t)./grid.g(gt(:,t+1));%b':=(1-cat*x')*bn/G'

        % Track CAT bond triggers
        cat_active = (cat(:,t) > 0);% CAT bond was bought
        triggered = cat_active & (xt(:,t+1) > dgrid.catthreshold);% CAT bond triggered

        cat_trigger_total = cat_trigger_total + cat_active;%count exposure
        cat_trigger_count = cat_trigger_count + triggered;%count actual triggers

    else
        bt(:,t+1)= bn(:,t)./grid.g(gt(:,t+1));
    end

    byt(:,t+1)=bt(:,t+1)./( exp(zt(:,t+1)).*kt(:,t+1).^par.alpha );%by':=b'/(k'^alpha)
    t;
    if byt(:,t+1)>par.costy*exp(zt(:,t+1)).*grid.g(gt(:,t+1)).^par.psi
        def(:,t+1)=1;
    end
    
    mt(:,t+1)=(1-def(:,t+1)*par.costy.*grid.g(gt(:,t+1)).^par.psi).*exp(zt(:,t+1)).*kt(:,t+1).^par.alpha - (1-def(t+1))*bt(:,t+1) + (1-par.delta)*kt(:,t+1);    
    
    ct(:,t)=mt(t)-kn(:,t)+qt(:,t).*bn(:,t);
    nx(:,t)=bt(t)-qt(:,t).*bn(:,t);
end
out.v=vt;
out.k=kt;
out.by=byt;
out.b=bt;
out.byn = byn;
out.bn = bn;
out.def=def;
out.s=st;
out.q=qt;
out.kn = kn;
out.m = mt;
out.nx=nx;
out.c=ct;
out.cat=cat;
out.cat_trigger_count = cat_trigger_count;
out.cat_trigger_total = cat_trigger_total;
out.cat_trigger_rate = cat_trigger_count ./ max(cat_trigger_total, 1);%Avoid divide-by-zero
end