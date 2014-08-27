library(stabledist)     # stable
library(sde)            # OU
library(sn)             # skew-normal

######
# normal-looking things
######

# don't expect anything here
# quick_sim(t=1,n=1000,theta=0.1)


# too many shifts erases kurtosis
# quick_sim(t=100,n=1000,theta=10)

# small theta, low rate

######
# non-normal looking things
######

# theta-alpha interactions
# quick_sim(t=100,n=1000,opt_kernel=rstable,alpha=1.8,gamma=1,beta=0,delta=0)

# prob from CPP-type effects of zero jumps (small rate/brlen)
# quick_sim(t=100,n=1000,theta=1) 


######
# variables
######

# x0           : starting value, X(0)=x0
# t            : ending time of X(t)
# n            : number of samples from X(t)
# rate         : Poisson rate
# sigma        : OU variance parameter
# theta        : OU mean-reversion parameter
# opt_kernel   : OU mean parameter, Y(t) ~ opt_kernel(...), t ~ Poisson(rate)
# opt_mean_fix : Draw optima from noise distribution or from stoch. proc.?
# ...          : Passed to environmental optimum kernel

quick_sim = function(x0=0, t=1, n=1000, rate=3, theta=3, sigma=sqrt(3), opt_mean_fix=TRUE, opt_kernel=rnorm, ...)
{
    xx=sim_ou(x0,t,n,rate,theta,sigma,opt_mean_fix,opt_kernel,...)
    yy=as_traits(xx)
    hist(yy,breaks=n/10)
}

# zero-mean rsn
rsn0 = function(n,delta,alpha,omega)
{
    delta = alpha/sqrt(1+alpha^2)
    xi0 = -(omega*delta*sqrt(2/pi))
    return(rsn)
}

sim_ou = function(x0=0, t=1, n=1000, rate=1, theta=1, sigma=1, opt_mean_fix=TRUE, opt_kernel=rnorm, ...)
{
    xt = list()
    for (i in 1:n)
    {
        xt[[i]] = list()
        xt[[i]]$times = c()
        xt[[i]]$vals = c()
        xt[[i]]$opts = c()

        # iteration variables
        t_cur = 0.0
        x_cur = x0
        o_cur = 0.0
        j = 1
        dx = 0
        dt = 0

        update = TRUE
        while (update)
        {
            # sample new time
            dt = rexp(n=1,rate=rate)
            if (t_cur+dt > t)
            {
                update = FALSE
                dt = t - t_cur
            }
            t_cur = t_cur + dt

            # sample new trait
            x_cur = rcOU(n=1,Dt=dt,x0=x_cur,theta=c(o_cur*theta,theta,sigma))

            # sample new optimum
            if (update)
            {
                if (opt_mean_fix)
                    o_cur = opt_kernel(n=1,...)
                else
                    o_cur = o_cur+opt_kernel(n=1,...)
            }

            # store new time, trait, and optimum
            xt[[i]]$times[j] = t_cur
            xt[[i]]$vals[j] = x_cur
            xt[[i]]$opts[j] = o_cur

            j = j + 1
        }
    }
    return(xt)
}

as_traits = function(x)
{
    y=c()
    for (i in 1:length(x))
    {
        n=length(x[[i]]$vals)
        v=x[[i]]$vals[n]
        y=c(y,v)
    }
    return(y)
}

opt_stable = function(x)
{

    # sas
    fn = function(par,...)
    {
        ret = -sum(dstable(x,alpha=par[1],beta=0,gamma=par[2],delta=0,log=T))
        return(ret)
    }
    theta = optim(par = c(2*runif(1),rexp(1)),
                   fn = fn,
                lower = c(1e-6,1e-6),
                upper = c(2,1e6))

    return(theta)
}

opt_norm = function(x)
{

    # norm
    fn = function(par,...)
    {
        ret = -sum(dnorm(x,mean=0,sd=par[1],log=T))
        return(ret)
    }
    theta = optim(par = c(rexp(1)),
                   fn = fn,
                lower = c(1e-6),
                upper = c(1e6))

    return(theta)
}

opt_like_ratio = function(x)
{
    theta1 = opt_stable(x)
    p1 = sum(dstable(x, alpha=theta1$par[1], gamma=theta1$par[2], beta=0, delta=0, log=T))

    theta2 = opt_norm(x)
    p2 = sum(dnorm(x, mean=0, sd=theta2$par[1], log=T))

    return(p1-p2)
}

# using sde's rcOU instead
rou = function(n=1000,t=1.0,mu=0.0,theta=0.1,sigma=1.0,x0=0.0)
{
    dw = rnorm(n,0,sqrt(t/n))
    dt = t/n
    x = c(x0)
    for (i in 2:(n+1))
        x[i] = x[i-1] + theta*(mu-x[i-1])*dt + sigma*dw[i-1]
    return(x)
}

