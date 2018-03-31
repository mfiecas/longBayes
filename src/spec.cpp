#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;

const int nobs=134;
const vec freq=linspace<vec>(0,1-1.0/nobs,nobs);
const vec w=linspace<vec>(0,1,nobs+1);
const int nhf=68;
const vec fold=join_cols(join_cols(ones(1),2*ones(nhf-2)),ones(1));

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat bases(int J=6)
{
    vec Jvec=linspace<vec>(1,J,J);
    mat tmp(nhf,J);
    for(uword k=0;k<nhf;k++)
    {
        for(uword j=0;j<J;j++) tmp(k,j)=freq(k)*Jvec(j);
    }
    return(sqrt(2.0)*cos(2.0*(datum::pi)*tmp));
}
// [[Rcpp::depends(RcppArmadillo)]]
vec powv(double x, vec p)
{
    int n=p.n_elem;
    vec res(n);
    for(uword i=0;i<n;i++) res(i)=pow(x,p(i));
    return(res);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec pg(vec ts)//nhf
{
    vec tmp=square( abs(fft(ts)) )/nobs;
    return(tmp.head(nhf));
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec lpg(vec ts)//nhf
{
    return(log( pg(ts) ));
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec nsp(vec spec, double TR=3.0)//nhf-nhf
{
    vec spec0=join_cols(spec,flipud(spec.head(nhf-1)));
    double fnorm=accu((w.tail(nobs)-w.head(nobs))%(spec0.tail(nobs)+spec0.head(nobs)))/2;
    return(TR*spec/fnorm);
}
// [[Rcpp::depends(RcppArmadillo)]]
vec ksm(vec pgram, int m)//nhf to nhf
{
    vec weight=ones(2*m+1)/2.0/m;
    weight(0)/=2.0;weight(2*m)/=2.0;
    
    vec hd=pgram.subvec(1,m);
    vec tl=pgram.subvec(nhf-1-m,nhf-2);
    vec extpgram=join_cols(flipud(hd),join_cols(pgram,flipud(tl)));
    
    vec res(nhf);
    for(uword i=0;i<nhf;i++) res(i)=dot(weight,extpgram.subvec(i,i+2*m));
    return(res);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec atsm(vec pgram)//nhf to nhf
{
    uvec nz=find(pgram);
    int nnz=nz.n_elem;
    
    vec weight=ones(nnz);
    weight(0)=0.5;weight(nnz-1)=0.5;
    
    double gval=1e10;
    vec ksres=pgram;
    vec ksnew(nhf);
    vec pg_sm(nhf);
    vec npg_sm(nnz);
    double gnew;
    
    for(int m=5;m<=15;m++)
    {
        ksnew=ksm(pgram,m);
        pg_sm=pgram/ksnew;
        npg_sm=pg_sm(nz);
        gnew=dot(weight,(-log(npg_sm)+npg_sm-1))/nnz/(1-1/2.0/m)/(1-1/2.0/m);
        if(gnew<gval)
        {
            gval=gnew;
            ksres=ksnew;
        }
    }
    return(ksres);
}
// [[Rcpp::depends(RcppArmadillo)]]
mat Reta(double rho, vec eta)
{
    int V=eta.n_elem+1;
    mat M=eye(V,V);
    for(uword i=2;i<=V;i++)
    {
        for(uword j=1;j<=(i-1);j++)
        {
            M(i-1,j-1)=pow(rho,sum(eta.subvec(j-1,i-2)));
            M(j-1,i-1)=M(i-1,j-1);
        }
    }
    return(M);
}
// [[Rcpp::depends(RcppArmadillo)]]
double detR(double rho, vec eta)
{
    return(prod(1-powv(rho,2*eta)));
}
// [[Rcpp::depends(RcppArmadillo)]]
mat Rinv(double rho, vec eta)//only when V>1
{
    int V=eta.n_elem+1;
    vec tmp=(1+powv(rho,2*eta))/(1-powv(rho,2*eta));
    vec mid=(tmp.head(V-2)+tmp.tail(V-2))/2;
    vec di=join_cols(join_cols(ones(1)/(1-pow(rho,2*eta(0))),mid),ones(1)/(1-pow(rho,2*eta(V-2))));
    vec di1=-powv(rho,eta)/(1-powv(rho,2*eta));
    return(diagmat(di)+diagmat(di1,1)+diagmat(di1,-1));
}
// [[Rcpp::depends(RcppArmadillo)]]
vec Rinvb(vec b, vec mu, double rho, vec eta, int J=6)//common term
{
    int V=b.n_elem/J;
    if(V==1) return(b-mu);
    
    mat bmat=reshape(b, J, V);
    for(uword v=0;v<V;v++) bmat.col(v)-=mu;
    mat gam=Rinv(rho,eta)*bmat.t();
    return(vectorise(gam.t()));
}
// [[Rcpp::depends(RcppArmadillo)]]
double bRinvb(vec b, vec mu, double rho, vec eta, int J=6)//common sum term
{
    int V=b.n_elem/J;
    if(V==1) return(dot(b-mu,b-mu));
    
    mat bmat=reshape(b, J, V);
    for(uword v=0;v<V;v++) bmat.col(v)-=mu;
    mat tbmat=bmat.t();
    return(dot(vectorise(tbmat),vectorise(Rinv(rho,eta)*tbmat)));
}

// [[Rcpp::depends(RcppArmadillo)]]
double lud(vec q, vec mu, double tau, double rho, vec eta, vec y, double siga=1000)
{
    int V=y.n_elem/nhf;
    int J=q.n_elem/V-1;
    mat bas=bases(J);
    //int V=q.n_elem/(J+1);
    vec a=q.head(V);
    vec b=q.tail(V*J);
    vec g(nhf);
    vec tmp(nhf);
    double py=0;
    for(uword v=1; v<=V; v++)
    {
        g=a(v-1)+bas*b.subvec(J*(v-1),J*v-1);
        tmp=g+exp(y.subvec(nhf*(v-1),nhf*v-1)-g);
        py+=dot(fold,tmp);
    }
    py/=(-2);
    double pa=-dot(a,a)/2/siga;
    double pb=-bRinvb(b,mu,rho,eta,J)/2/tau;
    return(py+pa+pb);
}
// [[Rcpp::depends(RcppArmadillo)]]
vec grad(vec q, vec mu, double tau, double rho, vec eta, vec y, double siga=1000)
{
    int V=y.n_elem/nhf;
    int J=q.n_elem/V-1;
    mat bas=bases(J);
    vec a=q.head(V);
    vec b=q.tail(V*J);
    vec g(nhf);
    vec tmp(nhf);
    vec pya(V);
    vec pyb(size(b));
    for(uword v=1; v<=V; v++)
    {
        g=a(v-1)+bas*b.subvec(J*(v-1),J*v-1);
        tmp=1-exp(y.subvec(nhf*(v-1),nhf*v-1)-g);
        pya(v-1)=-(dot(fold,tmp))/2;
        pyb.subvec(J*(v-1),J*v-1)=-bas.t()*(fold%tmp)/2;
    }
    vec pa=-a/siga;
    vec pb=-Rinvb(b,mu,rho,eta,J)/tau;
    return(join_cols(pya+pa,pyb+pb));
}
// [[Rcpp::depends(RcppArmadillo)]]
double ludrho(double rho, vec eta, vec b, vec mu, double tau, int J=6)
{
    return(-J*log(detR(rho,eta))/2-bRinvb(b,mu,rho,eta,J)/2/tau);
}
// [[Rcpp::depends(RcppArmadillo)]]
double gradrho(double rho, vec eta, vec b, vec mu, double tau, int J=6)
{
    int V=b.n_elem/J;
    if (V==1) return(0);
    
    vec r2e=powv(rho,2*eta);
    vec d1=eta%powv(rho,2*eta-1)/(1-r2e);
    vec d2=d1/(1-r2e);
    vec c=d2%(1+r2e)/powv(rho,eta);
    mat bmat=reshape(b, J, V);
    for(uword v=0;v<V;v++) bmat.col(v)-=mu;
    vec bno=sum(square(bmat)).t();
    vec bcr=sum(bmat.cols(0,V-2)%bmat.cols(1,V-1)).t();
    
    double val1=sum(d1);
    double val2=dot(d2,bno.head(V-1)+bno.tail(V-1));
    double val3=dot(c,bcr);
    return(J*val1-val2/tau+val3/tau);
}
// [[Rcpp::depends(RcppArmadillo)]]
vec abstep(vec q_cur, double epsilon, uword L, vec mu, double tau, double rho, vec eta, vec y, double siga=1000)
{
    vec q=q_cur;
    vec p=randn<vec>(size(q));
    vec p_cur=p;
    
    p+=epsilon*grad(q,mu,tau,rho,eta,y,siga)/2;
    for(uword i=1; i<=L; i++)
    {
        q+=epsilon*p;
        if (i!=L) p+=epsilon*grad(q,mu,tau,rho,eta,y,siga);
    }
    p+=epsilon*grad(q,mu,tau,rho,eta,y,siga)/2;
    
    p=-p;
    
    double lud_cur=lud(q_cur,mu,tau,rho,eta,y,siga);
    double K_cur=sum(square(p_cur))/2;
    double lud_pro=lud(q,mu,tau,rho,eta,y,siga);
    double K_pro=sum(square(p))/2;
    
    if (R::runif(0,1) < exp(-lud_cur+lud_pro+K_cur-K_pro)) return(q);
    else return(q_cur);
}
// [[Rcpp::depends(RcppArmadillo)]]
vec mustep(vec b, double tau, double rho, vec eta, int J=6, double sigmu=1000)
{
    int V=b.n_elem/J;
    double va;vec ct;
    if(V==1){
        va=1/(1/tau+1/sigmu);
        ct=b/tau*va;
    }
    else{
        mat bmat=reshape(b, J, V);
        vec sumRinv=sum(Rinv(rho,eta),1);
        va=1.0/(sum(sumRinv)/tau+1/sigmu);
        ct=bmat*sumRinv/tau*va;
    }
    return(ct+sqrt(va)*randn<vec>(J));
}
// [[Rcpp::depends(RcppArmadillo)]]
double taustep(vec b, vec mu, double rho, vec eta, double Ltau=1000, int J=6)
{
    double tau;
    do {
        tau=1/R::rgamma(0.5*(b.n_elem-1),2.0/bRinvb(b,mu,rho,eta,J));//default is scale
    } while (tau>Ltau);
    return(tau);
}
// [[Rcpp::depends(RcppArmadillo)]]
double rhostep(double rho_cur, vec eta, double epsilon, uword L, vec b, vec mu, double tau, int J=6)
{
    double rho=rho_cur;
    double p=R::rnorm(0,1);
    double p_cur=p;
    
    double rhotmp;
    double ptmp;
    
    p+=epsilon*gradrho(rho,eta,b,mu,tau,J)/2;
    for(uword i=1; i<=L; i++)
    {
        ptmp=p;
        rhotmp=rho+epsilon*ptmp;
        while( rhotmp<=0 || rhotmp>=1 )
        {
            ptmp=-ptmp;
            if(rhotmp>=1) rhotmp=2-rhotmp;
            else rhotmp=-rhotmp;
        }
        rho=rhotmp;
        p=ptmp;
        
        if (i!=L) p+=epsilon*gradrho(rho,eta,b,mu,tau,J);
    }
    p+=epsilon*gradrho(rho,eta,b,mu,tau,J)/2;
    
    p=-p;
    
    double lud_cur=ludrho(rho_cur,eta,b,mu,tau,J);
    double K_cur=p_cur*p_cur/2;
    double lud_pro=ludrho(rho,eta,b,mu,tau,J);
    double K_pro=p*p/2;
    
    if (R::runif(0,1) < exp(-lud_cur+lud_pro+K_cur-K_pro)) return(rho);
    else return(rho_cur);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat mcmc(vec y, vec eta, vec q, double leps=0.065, uword L=15, double lepsrho=0.13, uword Lrho=10, uword niter=2e4, double Ltau=1000, double siga=1000, double sigmu=1000)
{
    int V=y.n_elem/nhf;
    int J=q.n_elem/V-1;
    double tau=10;
    double rho=0.5;
    vec mu=zeros<vec>(J);
    q=abstep(q,leps*R::runif(0,1),L,mu,tau,rho,eta,y,siga);
    vec b=q.tail(V*J);
    mu=mustep(b,tau,rho,eta,J,sigmu);
    tau=taustep(b,mu,rho,eta,Ltau,J);
    if(V>1) rho=rhostep(rho,eta,lepsrho*R::runif(0,1),Lrho,b,mu,tau,J);
    
    mat q_out=q.t();
    mat mu_out=mu.t();
    NumericVector tau_out(1); tau_out(0)=tau;
    NumericVector rho_out(1); rho_out(0)=rho;
    
    for(uword i=1; i<niter; i++)
    {
        q=abstep(q,leps*R::runif(0,1),L,mu,tau,rho,eta,y,siga);
        vec b=q.tail(V*J);
        mu=mustep(b,tau,rho,eta,J,sigmu);
        tau=taustep(b,mu,rho,eta,Ltau,J);
        if(V>1) rho=rhostep(rho,eta,lepsrho*R::runif(0,1),Lrho,b,mu,tau,J);
        
        q_out.insert_rows(i, q.t());
        mu_out.insert_rows(i, mu.t());
        tau_out.push_back(tau);
        if(V>1) rho_out.push_back(rho);
    }
    if(V>1) return( join_rows(q_out,join_rows(mu_out,join_rows(as<vec>(tau_out),as<vec>(rho_out))) ) );
    else return( join_rows(q_out, join_rows(mu_out,join_rows(as<vec>(tau_out), zeros(niter))) ) );
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat q_lspec(vec q, int J=6)
{
    mat bas=bases(J);
    int V=q.n_elem/(J+1);
    vec a=q.head(V);
    vec b=q.tail(V*J);
    vec g(nhf);
    mat out(nhf,V);
    for(uword v=1; v<=V; v++)
    {
        g=a(v-1)+bas*b.subvec(J*(v-1),J*v-1);
        out.col(v-1)=g;
    }
    return(out);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec mc_q(mat mcsamp, int J=6, int burn=0)//half spec
{
    int nc=mcsamp.n_cols;
    mcsamp.shed_cols(nc-J-2,nc-1);
    if (burn>0) mcsamp.shed_rows(0,burn-1);
    vec qmean=mean(mcsamp).t();
    return(qmean);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat lspec(mat mcsamp, int J=6, int burn=0)//nhf
{
    vec qmean=mc_q(mcsamp,J,burn);
    return(q_lspec(qmean,J));
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat tsl_mc(mat tsl, int J=6)//tsl has age in first row
{
    int V=tsl.n_cols;
    vec eta;
    if(V==1) eta=zeros(1);
    else
    {
        vec age=trans(tsl.row(0));
        eta=4*(age.tail(V-1)-age.head(V-1));//time unit is 3 months
    }
    tsl.shed_row(0);
    
    mat yl(nhf,V);
    for(uword v=1; v<=V; v++)
    {
        yl.col(v-1)=lpg(tsl.col(v-1));
    }
    vec y=vectorise(yl);
    mat mcsamp=mcmc(y,eta,0.5*ones<vec>(V*(J+1)));
    return(mcsamp);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec tsl_q(mat tsl, int J=6, int burn=0)//tsl has age in first row
{
    mat mcsamp=tsl_mc(tsl,J);
    vec qmean=mc_q(mcsamp,J,burn);
    return(qmean);
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat lspl(mat tsl, int J=6, int burn=0)//tsl has age in first row
{
    mat mcsamp=tsl_mc(tsl,J);
    return(lspec(mcsamp,J,burn));
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat spl(mat tsl, int J=6, int burn=0)//tsl has age as first row
{
    return(exp( lspl(tsl,J,burn)) );
}
