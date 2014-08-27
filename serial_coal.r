source("jumpdif.r")

library(ape)
library(ggplot2)
library(reshape2)
library(scales)
library(parallel)
library(sn)
library(VGAM)
library(moments)
library(diptest)
library(parallel)
require(gridExtra)
library(gtable)

###################
#  simulate data  #
###################
sim_qts = function(nrep=50,loc=1:8,sam=1:8,npop=2000,theta=.1,kernel=rnorm,scale=1,alpha=2.0,alpha_sn=0.0,xi_sn=0.0,ncores=16)
{
    # loci
    b = list()
    nloc=length(loc)
    nsam=length(sam)
    loc_offset=-1
    for (n in loc)
    {
        n_ = n+loc_offset
        b[[n]] = list()
        # samples
        for (m in sam)
        {
            scale_per_locus = scale/((2^n_)^(1/alpha))
            if (identical(kernel,rnorm))
            {
                print(c(n,m))
                b[[n]][[m]] = batch_ms(nrep=nrep,spop=2^m,npop=npop,nloci=2^n_,ncores=ncores,theta=theta,kernel=kernel,sd=scale_per_locus)
            }
            else if (identical(kernel,rsn))
            {
                delta_sn = alpha_sn/sqrt(1+alpha_sn^2)
                omega_sn = scale_per_locus / sqrt(1-(2*delta_sn^2)/pi)
                xi_0_sn = -(omega_sn * delta_sn * sqrt(2/pi))
                b[[n]][[m]] = batch_ms(nrep=nrep,spop=2^m,npop=npop,nloci=2^n_,ncores=ncores,theta=theta,kernel=kernel,xi=xi_0_sn+xi_sn,omega=omega_sn,alpha=alpha_sn)
            }
            else if (identical(kernel,rlaplace))
            {
                print(c(n,m))
                b[[n]][[m]] = batch_ms(nrep=nrep,spop=2^m,npop=npop,nloci=2^n_,ncores=ncores,theta=theta,kernel=kernel,location=0,scale=scale_per_locus/sqrt(2))
            }
            else if (identical(kernel,rstable))
            {
                b[[n]][[m]] = batch_ms(nrep=nrep,spop=2^m,npop=npop,nloci=2^n_,ncores=ncores,theta=theta,kernel=kernel,gamma=scale_per_locus,alpha=alpha,beta=0)
            }
            else
            {
                print("kernel must equal \"rnorm\", \"rsn\", \"rlaplace\", or \"rstable\".")
                return
            }

        }
    }
    kernel_fn = c("rnorm","rsn","rstable","rlaplace")
    kernel_str = c("normal","skew-normal","stable","Laplace")
    kernel_query=as.character(substitute(kernel))
    kidx = which(kernel_fn==kernel_query)
    ks = kernel_str[kidx]
    if (kidx==2)
        ks = paste(ks,"-",solve_alpha_for_skew(alpha_sn),sep="")
    else if (kidx==3)
        ks = paste(ks,"-",alpha,sep="")
    b$kernel_str=ks
    b$kernel_fn=kernel
    b$skew=solve_delta_for_skew(solve_alpha_for_delta(alpha_sn))

    return(b)
}

# create ms command string
make_cmd_str = function( out_file="ms_trees.txt", ms_exec="~/apps/msdir/ms",spop=2,nloci=10, theta=10.0,div_time=0.1,seeds=sample(1:10000,3))
{
    cmd_str = ms_exec
    #paste(cmd_str,sum(spop),nloci,"-T","-t",theta,"-I",length(spop),spop[1],spop[2],0,"-ej",div_time/2,1,2,sep=" ")
    #paste(cmd_str,spop,nloci,"-t",theta,">",paste(getwd(),"/",out_file,sep=""),sep=" ")
    paste(cmd_str,spop,nloci,"-t",theta,"-seeds",seeds[1],seeds[2],seeds[3],sep=" ")
}

# runs batch jobs of ms to simulate data
batch_ms = function( out_file="ms_trees.txt",ms_path="~/apps/msdir/ms",nrep=1,spop=2,npop=2000,nloci=10,ncores=8,theta=10,div_time=0.1,kernel=rnorm,...)
{
    # redirect ms output to out_file
    simlist = list()
    dt = proc.time()[3]

    simlist = mclapply(1:nrep, function(josh) {
    #for (josh in 1:nrep) {
        t = dt
        dt = proc.time()[3]

        # read in ms output
        #ms_out = system(cmd_str, intern=TRUE)
        cmd_str = make_cmd_str(out_file,ms_path,npop,nloci,theta);

        print(cmd_str)
        ms_out = system(cmd_str,intern=T)
        ds = proc.time()[3]
        write(sprintf("Replicate %d, nl=%d, spop=%d, t=%f",josh,nloci,spop,as.numeric(t)),"")
        nrow = length(ms_out)
        sim = c()
        ss_mtx = list()
        locus_idx = 0
        rand_idx_list = sample(1:npop, spop, replace=F)
        sample_idx = 0
        subsample_idx = 0
        skip_lines = TRUE
        nsegsites = c()

        for (i in 1:nrow)
        {
            line = ms_out[i]
            first_char = substr(line,1,1)

            # new sim
            if (substr(line,1,2) == "//")
            {
                skip_lines = FALSE
                sample_idx = 0
                subsample_idx = 0
                locus_idx = locus_idx + 1
                #sim[[locus_idx]] = list()
            }
            # skip ms header
            else if (skip_lines == TRUE)
            {
                next
            }
            # tree string
            else if (first_char == "(")
            {
                # do nothing
            }
            # num segsites
            else if (first_char == "s")
            {
                nss = as.numeric(unlist(strsplit(line," "))[2])
                nsegsites[locus_idx] = nss
                ss_mtx[[locus_idx]] = matrix(ncol=nss,nrow=spop)
                sample_idx = 0
                subsample_idx = 0
            }
            # segsites
            else if (first_char == "0" || first_char == "1")
            {
                sample_idx = sample_idx + 1
                if ( !(sample_idx %in% rand_idx_list) )
                {
                    next
                }
                subsample_idx = subsample_idx + 1

                sample_segsites = as.numeric(unlist(strsplit(line,"")))
                ss_mtx[[locus_idx]][subsample_idx,] = sample_segsites
            }
        }

        # simulate against sampled alleles per locus
        sim$trait = numeric(spop)
        for (i in 1:nloci)
        {
            if (nsegsites[i] != 0)
            {
                mut_effects = kernel(n=nsegsites[i],...)
                sim$trait = sim$trait + ss_mtx[[i]]%*%mut_effects
            }
        }
        rm(ss_mtx)
        sim$ss_mtx = NULL
        sim$median_trait = sim$trait-median(sim$trait)
        sim$mean_trait = sim$trait-mean(sim$trait)
        #gc()
        return(sim)
    },mc.cores=ncores)
    return(simlist)
}

append_reps = function(...,loc=1:8,sam=10:10)
{
    sim_list = list(...)
    ret_sim = sim_list[[1]]
    nsim=length(sim_list)
    curr_idx = 1
    end_idx=length(ret_sim[[loc[1]]][[sam[1]]])

    for (i in 2:nsim)
    {
        curr_idx = curr_idx + end_idx
        end_idx = length(sim_list[[i]][[loc[1]]][[sam[1]]])
        for (j in 1:end_idx)
        {
            for (k in loc)
            {
                for (l in sam)
                {
                    ret_sim[[k]][[l]][[curr_idx+j]] = sim_list[[i]][[k]][[l]][[j]]
                }
            }
        }
    }
    return(ret_sim)
}


############
# KS tests #
############

compute_ks = function(b,loc=1:8,sam=1:8,nrep=50,theta=0.1,scale=1.0,alpha=2.0,p_val=0.05)
{
    nloc=length(loc)
    nsam=length(sam)
    ks_stat=array(0,dim=c(nloc,nsam,nrep),dimnames=c("nloc","nsam","nrep"))
    ks_pvalue=array(0,dim=c(nloc,nsam,nrep),dimnames=c("nloc","nsam","nrep"))
    ks_pvalue_sig = array(0,dim=c(nloc,nsam),dimnames=c("nloc","nsam"))
    for (n in loc)
    {
        for (m in sam)
        {
            n_p_sig = 0
            for (l in 1:nrep)
            {
                cur_traits = b[[n]][[m]][[l]]$trait
                med_traits = cur_traits-median(cur_traits)
                #print(cur_traits)
                # get KS from traits
                if (identical(b$kernel_fn,rnorm))
                    ks_res = ks.test(med_traits,pnorm,sd=sqrt(1/2*theta*scale^2))
                else if (identical(b$kernel_fn,rsn))
                    ks_res = ks.test(med_traits,pnorm,sd=sqrt(1/2*theta*scale^2))
                else if (identical(b$kernel_fn,rlaplace))
                    ks_res = ks.test(med_traits,pnorm,sd=sqrt(1/2*theta*scale^2))
                else if (identical(b$kernel_fn,rstable))
                    ks_res = ks.test(med_traits,pstable,gamma=scale*theta/2,alpha=alpha,beta=0)
                #print(ks_res)
                b[[n]][[m]][[l]]$ks_stat = ks_res$statistic
                b[[n]][[m]][[l]]$ks_pvalue = ks_res$p.value

                ks_stat[n,m,l] = b[[n]][[m]][[l]]$ks_stat
                ks_pvalue[n,m,l] = b[[n]][[m]][[l]]$ks_pvalue

                if (ks_res$p.value < p_val)
                    n_p_sig = n_p_sig + 1
            }
            ks_pvalue_sig[n,m] = n_p_sig / nrep
        }
    }
    b$ks_stat = ks_stat
    b$ks_pvalue = ks_pvalue
    b$ks_stat_mean = apply(ks_stat,1:2,mean)
    b$ks_pvalue_mean = apply(ks_pvalue,1:2,mean)
    b$ks_pvalue_freq = ks_pvalue_sig
    return(b)
}

plot_ksfp = function(b,title=b$kernel_str)
{
    dev.new()
    nloc=dim(b$ks_pvalue_freq)[1]
    nsam=dim(b$ks_pvalue_freq)[2]
    ggplot(melt(b$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title=title) +
        xlab("# loci") +
        ylab("# samples") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+1 )) +
        scale_fill_gradientn(name="freq(p<0.05)",limits=c(0.0,1.0),colours=c("white","black"),breaks=c(0.0,1.0)) +
        theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.y = element_text(size=16), axis.title.x = element_text(size=16))
}
plot_kspm = function(b,title=b$kernel_str)
{
    dev.new()
    ggplot(melt(b$ks_pvalue_mean),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title=title) +
        xlab("# loci") +
        ylab("# samples") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_fill_gradientn(name="1-p",limits=c(0.0,1.0),colours=c("white","black"),breaks=c(0.0,1.0)) +
        theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.y = element_text(size=16), axis.title.x = element_text(size=16))
}

plot_kssm = function(b, title=b$kernel_str)
{
    dev.new();
    ggplot(melt(0.5-b$ks_stat_mean),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        #geom_text(colour="white",aes(label=round(value,2))) +
        labs(title=title) +
        xlab("# loci") +
        ylab("# samples") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_fill_gradientn(name="sampling\ndistribution",limits=c(0.0,0.5),colours=c("white","black"),breaks=c(0.0,0.5),labels=format(c("non-\nnormal","normal"))) +
        theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.y = element_text(size=16), axis.title.x = element_text(size=16))
}

plot_all_ksfp_diff = function(k1,k2,k3)
{

    dev.new()
    nloc=dim(k1$ks_pvalue_freq)[1]
    nsam=dim(k1$ks_pvalue_freq)[2]

    p1 = ggplot(melt(k1$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="normal minus zero") +
        xlab("L") +
        ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+2 )) +
        scale_fill_gradientn(limits=c(-1.0,1.0),colours=c("red","white","black"))
        #theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())
        #theme(axis.text.y=element_blank(), axis.title.y=element_blank())

    p2 = ggplot(melt(k2$ks_pvalue_freq[,3:nsam]-k1$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="skew-normal (0.9) minus normal") +
        xlab("L") +
        #ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        #scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+1 )) +
        scale_fill_gradientn(limits=c(-1.0,1.0),colours=c("red","white","black")) +
        #theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())
        theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())


    p3 = ggplot(melt(k3$ks_pvalue_freq[,3:nsam]-k1$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="Laplace minus normal") +
        xlab("L") +
        #ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        #scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+1 )) +
        scale_fill_gradientn(name="freq(p<0.05)\nA minus B",limits=c(-1.0,1.0),colours=c("red","white","black"),breaks=c(-1.0,0.0,1.0)) +
        theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())

    legend = gtable_filter(ggplot_gtable(ggplot_build(p3)),"guide-box")

    pp = grid.arrange(arrangeGrob(  p1 + theme(legend.position="none"),
                                    p2 + theme(legend.position="none"),
                                    p3 + theme(legend.position="none"),
                                    nrow = 1),
                                legend,
                                widths=c(6,1),
                                nrow = 1)

    return(pp)

}

plot_all_ksfp = function(k1,k2,k3)
{

    dev.new()
    nloc=dim(k1$ks_pvalue_freq)[1]
    nsam=dim(k1$ks_pvalue_freq)[2]

    p1 = ggplot(melt(k1$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="normal") +
        xlab("L") +
        ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+2 )) +
        scale_fill_gradientn(limits=c(0.0,1.0),colours=c("white","black"))
        #theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank())
        #theme(axis.text.y=element_blank(), axis.title.y=element_blank())

    p2 = ggplot(melt(k2$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="skew-normal (0.9)") +
        xlab("L") +
        #ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        #scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+1 )) +
        scale_fill_gradientn(limits=c(0.0,1.0),colours=c("white","black")) +
        theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())


    p3 = ggplot(melt(k3$ks_pvalue_freq[,3:nsam]),aes(Var1,Var2,fill=value)) +
        geom_raster() +
        labs(title="Laplace") +
        xlab("L") +
        #ylab("n") +
        scale_x_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x-1 )) +
        #scale_y_continuous(breaks=c(seq(0:10)),labels=math_format(expr=2^.x,format=function(x) x+1 )) +
        scale_fill_gradientn(name="freq(p<0.05)",limits=c(0.0,1.0),colours=c("white","black"),breaks=c(0.0,1.0)) +
        theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())

    legend = gtable_filter(ggplot_gtable(ggplot_build(p3)),"guide-box")

    pp = grid.arrange(arrangeGrob(  p1 + theme(legend.position="none"),
                                    p2 + theme(legend.position="none"),
                                    p3 + theme(legend.position="none"),
                                    nrow = 1),
                                legend,
                                widths=c(6,1),
                                nrow = 1)

    return(pp)

}
###################
# compute moments #
###################

get_moment = function(b, loc=1:7, sam=10:10, moment=4, use_moments=TRUE, central=TRUE)
{
    vals = c()

    olveyLoc=TRUE
    if (length(loc)==1 && length(sam)>1)
        byLoc=FALSE

    nreps=sum(names(b[[loc[2]]][[sam[1]]])=="")
    if (nreps==0)
        nreps=length(b[[loc[2]]][[sam[1]]])
    reps = 1:nreps
    for (i in loc) {
        for (j in sam) {
            tmp = c()
            for (k in reps) {
                if (use_moments)
                {
                    tmp[k] = moment(b[[i]][[j]][[k]]$trait,order=moment,central=central)
                }
                else if (moment == 2 && !use_moments) {
                    tmp[k] = var(b[[i]][[j]][[k]]$trait)
                }

                else if (moment == 3 && !use_moments) {
                    tmp[k] = skewness(b[[i]][[j]][[k]]$trait,method="moment")
                }
                else if (moment == 4 && !use_moments) {

                    #tmp[k] = mykurtosis(b[[i]][[j]][[k]]$trait)
                    tmp[k] = kurtosis(b[[i]][[j]][[k]]$trait,method="moment")
                }
            }
            if (byLoc)
                vals[i] = mean(tmp,na.rm=T)
            else
                vals[j] = mean(tmp,na.rm=T)
        }
    }
    return(vals)
}

mykurtosis = function(x)
{
    n = length(x)
    x_bar = mean(x)
    k2 = sum( (x-x_bar)^2 ) / (n-1)
    a = ( (n+1)*n ) / ( (n-1)*(n-2)*(n-3) )
    b = sum( (x-x_bar)^4 ) / k2^2
    c = -3 * (n-1)^2 / ( (n-2)*(n-3) )
    return(a*b+c)
}

solve_skew_for_delta = function(skew,sign=1)
{
    if (skew>1) sign*1
    else sign*sqrt( (pi/2) * abs(skew)^(2/3) / ( abs(skew)^(2/3) + ((4-pi)/2)^(2/3) ) )
}

solve_delta_for_alpha= function(delta,sign=1)
{
    sign*delta/sqrt(1-delta^2)
}

solve_skew_for_alpha = function(skew, sign=1)
{
    solve_delta_for_alpha(solve_skew_for_delta(skew, sign),sign)
}

solve_alpha_for_delta=function(alpha)
{
    alpha/sqrt(1+alpha^2)
}

solve_delta_for_skew = function(delta)
{
    (4-pi)/2 * (delta*sqrt(2/pi))^3 / (1-2*delta^2/pi)^(3/2)
}

solve_alpha_for_skew=function(alpha)
{
    solve_alpha_for_delta(solve_delta_for_skew(alpha))
}

sn_skew = function(omega=1,delta=1.0)
{
    s = (4-pi)/2 * (delta*sqrt(2/pi))^3 / (1-2*delta^2/pi)^(3/2)
    v = omega^2 * (1-2*delta^2/pi)
    return( s*v^(3/2) )
}

sn_kurtosis = function(omega=1,delta=1) {
    k = 2*(pi-3)*(delta*sqrt(2/pi))^4/(1-2*delta^2/pi)^2
    v = omega^2*(1-2*delta^2/pi)
    (k+3)*v^2
}

#in all these functions, sigma is the TARGET sigma
second_moment_kernel = function(theta,sigma,...) {
    1/2*theta*sigma^2
}

third_moment_kernel = function(theta,sigma,L,kernel,skew=0) {
    if (identical(kernel,rnorm)) {
        mu3 = 0
    } else if (identical(kernel,rlaplace)) {
        mu3 = 0
    } else if (identical(kernel,rsn)) {
        delta = solve_skew_for_delta(skew)
        tau = sigma/sqrt(L)
        omega = tau/sqrt(1-2*delta^2/pi)
        mu3 = sn_skew(omega=omega,delta=delta)
    }
    1/6*theta*L*mu3
}

fourth_moment_kernel = function(theta,sigma,L,kernel,skew=0) {
    if (identical(kernel,rnorm)) {
        tau = sigma/sqrt(L)
        mu2 = tau^2
        mu4 = 3*tau^4
    } else if (identical(kernel,rlaplace)) {
        tau = sigma/sqrt(L)
        mu2 = tau^2
        mu4 = 6*tau^4
    } else if (identical(kernel,rsn)) {
        delta = solve_skew_for_delta(skew)
        tau = sigma/sqrt(L)
        omega = tau/sqrt(1-2*delta^2/pi)
        mu2 = omega^2*(1-2*delta^2/pi)
        mu4 = sn_kurtosis(omega,delta)
    }
    3/4*L^2*theta^2*mu2^2+1/4*L*theta*(2*theta*mu2^2+mu4)
}

get_expected_moment = function(kernel, loc=1:8, sam=10:10, theta=2, sigma=1, skew=0, moment=4)
{
    moment_fn=c(NA, second_moment_kernel, third_moment_kernel, fourth_moment_kernel)[[moment]]
    v=c()
    #nloc=2^(max(loc)-1)
    nloc=10*(max(loc))
    locseq=seq(0,max(loc)-1,by=0.1)
    for (i in 1:length(locseq)) v[i]=moment_fn(theta,sigma,2^(locseq[i]),kernel,skew)
    #for (i in loc) v[i]=moment_fn(theta,sigma,2^(i-1),kernel,skew)
    return(v)
}

combine_moments = function(...,loc=1:8,sam=10:10,theta=2,sigma=1,skew=0,moment=4)
{
    input_list = list(...)
    nloc=length(loc)
    nargs=length(input_list)

    moment_vals=data.frame()
    for (k in 1:nargs)
    {
        write(sprintf("Getting simulated moment: %s (%d of %d)",input_list[[k]]$kernel_str, k, nargs),"")
        v = get_moment(b=input_list[[k]],loc=loc,sam=sam,moment=moment)
        for (i in loc)
        {
            moment_vals[nloc*(k-1)+i,1]=paste("simulated_",input_list[[k]]$kernel_str,sep="")
            moment_vals[nloc*(k-1)+i,2]=i-1
            moment_vals[nloc*(k-1)+i,3]=v[i]
            moment_vals[nloc*(k-1)+i,4]="simulated"
            moment_vals[nloc*(k-1)+i,5]=input_list[[k]]$kernel_str
        }
    }

    start_row=dim(moment_vals)[1] # nloc*nargs ??
    for (k in 1:nargs)
    {
        write(sprintf("Getting analytical moment: %s (%d of %d)",input_list[[k]]$kernel_str, k, nargs),"")
        ev = get_expected_moment(kernel=input_list[[k]]$kernel_fn,loc=loc,sam=sam,theta=theta,sigma=sigma,moment=moment,skew=input_list[[k]]$skew)
        nev=length(ev)
        for (i in 1:length(ev))
        {
            moment_vals[start_row+nev*(k-1)+i,1]=paste("analytical_",input_list[[k]]$kernel_str,sep="")
            #moment_vals[start_row+nev*(k-1)+i,2]=i
            moment_vals[start_row+nev*(k-1)+i,2]=(i-1)*0.1
            moment_vals[start_row+nev*(k-1)+i,3]=ev[i]
            moment_vals[start_row+nev*(k-1)+i,4]="analytical"
            moment_vals[start_row+nev*(k-1)+i,5]=input_list[[k]]$kernel_str
        }
    }
    return(moment_vals)
}

plot_moments = function(df,moment=4,legend=TRUE)
{
    title_str=c("1st moment", "2nd moment", "3rd moment", "4th moment")
    #dev.new()
    p = ggplot(df, aes(x=V2,y=V3,color=V5,group=V1,linetype=V4)) +
        geom_line() +
        xlab("L") +
        scale_linetype_manual(breaks=c("simulated","analytical"),values=c(2,1)) +
        scale_x_continuous(breaks=c(seq(0,length(unique(df$V2)))),labels=math_format(2^.x)) +
        guides(color=guide_legend(title="kernel"), linetype=guide_legend(title="type"))
        #theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title.y = element_text(size=16), axis.title.x = element_text(size=16))

    #if (legend) p = p + guides(color=guide_legend(title="kernel"), linetype=guide_legend(title="type"))
    #else p = p + theme(legend.position="none")

    if (moment==2) p = p + ylab(expression(h[2]))
    if (moment==3) p = p + ylab(expression(h[3]))
    if (moment==4) p = p + ylab(expression(h[4]))

    return(p)
}


plot_all_moments = function(m2,m3,m4)
{

    dev.new()
    p2 = ggplot(m2, aes(x=V2,y=V3,color=V5,group=V1,linetype=V4)) +
        geom_line() +
        scale_linetype_manual(breaks=c("simulated","analytical"),values=c(2,1)) +
        scale_x_continuous(breaks=c(seq(0,length(unique(m4$V2)))),labels=math_format(2^.x)) +
        #theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank()) +
        guides(color=guide_legend(title="kernel"), linetype=guide_legend(title="type")) +
        ylab(expression(h[2])) +
        xlab("L")

    p3 = ggplot(m3, aes(x=V2,y=V3,color=V5,group=V1,linetype=V4)) +
        geom_line() +
        scale_linetype_manual(breaks=c("simulated","analytical"),values=c(2,1)) +
        scale_x_continuous(breaks=c(seq(0,length(unique(m4$V2)))),labels=math_format(2^.x)) +
        #theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank()) +
        guides(color=guide_legend(title="kernel"), linetype=guide_legend(title="type")) +
        ylab(expression(h[3])) +
        xlab("L")

    p4 = ggplot(m4, aes(x=V2,y=V3,color=V5,group=V1,linetype=V4)) +
        geom_line() +
        scale_linetype_manual(breaks=c("simulated","analytical"),values=c(2,1)) +
        scale_x_continuous(breaks=c(seq(0,length(unique(m4$V2)))),labels=math_format(2^.x)) +
        guides(color=guide_legend(title="kernel"), linetype=guide_legend(title="type")) +
        ylab(expression(h[4])) +
        xlab("L")

    legend = gtable_filter(ggplot_gtable(ggplot_build(p2)),"guide-box")

    pp = grid.arrange(arrangeGrob(  p2 + theme(legend.position="none"),
                                    p3 + theme(legend.position="none"),
                                    p4 + theme(legend.position="none"),
                                    nrow = 1),
                                legend,
                                widths=c(6,1),
                                nrow = 1)

    return(pp)

}

g_legend = function(a.gplot) {
    tmp=ggplot_gtable(ggplot_build(a.gplot))
    leg=which(sapply(tmp$grobs, function(x) x$name)=="guide-box")
    legend=tmp$grobs[[leg]]
    return(legend)
}


######################
# test multimodality #
######################

# stores no. of replicates that are sig. non-unimodal under dip test
batch_diptest = function(b,loc=1:8,sam=10:10,nrep=50,p=0.05,d_th=0.05)
{
    # Bonferonni
    p_corr=p#/(length(loc)*length(sam)*nrep)

    for (i in loc)
    {
        for (j in sam)
        {
            pv=c()
            dv=c()
            n=0
            for (k in seq(nrep))
            {
                dt=dip.test(b[[i]][[j]][[k]]$trait)
                pv[k]=dt$p.value
                dv[k]=dt$statistic
                if (pv[k]<p_corr)# && dv[k]>d_th)
                    n=n+1
            }
            b[[i]][[j]]$f_dip = n/nrep
            b[[i]][[j]]$n_dip = n
            b[[i]][[j]]$dip = cbind(dv,pv)
        }
    }
    return(b)
}

combine_diptest = function(...,loc=1:8,sam=10:10)
{
    input_list = list(...)
    nloc=length(loc)
    nargs=length(input_list)

    dipvals=data.frame()
    for (k in 1:nargs)
    {
        write(sprintf("Combining %s (%d of %d)",input_list[[k]]$kernel_str, k, nargs))
        for (i in loc)
        {
            for (j in sam)
            {
                dipvals[nloc*(k-1)+i,1]=input_list[[k]]$kernel_str
                dipvals[nloc*(k-1)+i,2]=i-1
                dipvals[nloc*(k-1)+i,3]=input_list[[k]][[i]][[j]]$f_dip
            }
        }
    }

    return(dipvals)
}

plot_diptest = function(dfs)
{
    dev.new()
    ggplot(dfs, aes(x=V2,y=V3,color=V1,group=V1)) +
        geom_line() +
        scale_x_continuous(breaks=c(seq(0,length(unique(dfs$V2)))),labels=math_format(2^.x)) +
        xlab("L") +
        ylab("freq. multimodal") +
        guides(color=guide_legend(title="kernel"))
}



##################
##################
# CODE GRAVEYARD #
##################
##################

################################
# optimize/model fitting (old) #
################################

optim_coal = function(X,dkern=dnorm)
{
    par=c() # init parameter value
    lower=c() # lower bound
    upper=c() # upper bound

    # Josh probably has a better way to do this using elipsis magic
    if (identical(dkern,dnorm))
    {
        par=c(rexp(1))
        lower=c(10^-6)
        upper=c(10^6)
        fn = function(par,...)
        {
            ret = -sum(dnorm(X,mean=0,sd=par[1],log=T))
            return(ret)
        }
    }
    else if (identical(dkern,dstable))
    {
        par=c(2*runif(1),rexp(1))
        lower=c(0,0)
        upper=c(2,10^6)
        fn = function(par,...)
        {
            ret = -sum(dstable(X,alpha=par[1],beta=0,gamma=par[2],delta=0,log=T))
            return(ret)
        }
    }

    return(optim(   par=par,
                    fn=fn,
                    lower=lower,
                    upper=upper
            ))
}

# Can't remember if this is exactly what we wanted:
# Provided a kernel and data, find outliers wrt MLE
outlier_coal = function(X,p=5*10^-2,dkern=dnorm)
{
    optim_obj=optim_coal(X,dkern=dnorm)
    X_cdf = c()
    neg_tail_idx = c()
    pos_tail_idx = c()

    # get cdf from kernel
    if (identical(dkern,dnorm))
    {
        sd=optim_obj$par[1]
        X_cdf=pnorm(X,mean=0,sd=sd)
    }
    else if (identical(dkern,dstable))
    {
        alpha=optim_obj$par[1]
        gamma=optim_obj$par[2]
        X_cdf=pstable(X,alpha=alpha,gamma=gamma,beta=0,delta=0)
    }
    neg_tail_idx=which(X_cdf < p/2)
    pos_tail_idx=which(X_cdf > 1.0 - p/2)

    ret = list()
    ret$X = X
    ret$p = p
    ret$optim_obj = optim_obj
    ret$X_cdf = X_cdf
    ret$neg_tail_idx = neg_tail_idx
    ret$pos_tail_idx = pos_tail_idx

    # e.g. gives empirical threshhold value for left tail, approaching from -Inf
    # max(ret$X[ret$neg_tail_idx])

    return(ret)
}

# optim() e.g.
# foo.unconstr <- function(par, x) -sum(dnorm(x, par[1], par[2], log=TRUE))
# x <- rnorm(100,1,1)
# optim(c(1,1), foo.unconstr, lower=c(0,0), upper=c(5,5), method="L-BFGS-B", x=x)

int_ecf_dist = function(ecf,tcf,wtf,...)
{
    my_fn = ecf_dist(ecf,tcf,wtf,...)
    v=integrate(my_fn,-Inf,Inf)
    return(v)
}

int_ecf_dist_scale = function(ecf,tcf,wtf,gamma,...)
{
    my_fn = ecf_dist_scale(ecf,tcf,wtf,gamma,...)
    v=integrate(my_fn,-Inf,Inf)
    return(v)
}

ecf_dist_scale = function(ecf,tcf,wtf,gamma,...)
{
    #tcf should have gamma = 1, so it is standard
    return ( function(k) { (Mod(ecf(k,gamma) - tcf(k,...)))^2 * wtf(k)  } )
}

ecf_dist = function(ecf,tcf,wtf,...)
{
    return ( function(k) { (Mod(ecf(k) - tcf(k,...)))^2 * wtf(k)  } )
}

ecf_fn = function(x)
{
    return ( function(k) { rowMeans(exp(complex(imaginary=1) * k %*% t(x) )) } )
}

ecf_scale_fn = function(x)
{
    return ( function(k,gamma) { rowMeans(exp(complex(imaginary=1) * k %*% (t(x)/gamma) )) } )
}

cstable = function(k,alpha=2.0,beta=0.0,delta=0.0,gamma=1.0)
{
    im = complex(imaginary=1)
    if (alpha != 1)
        return(exp( -gamma^alpha * abs(k)^alpha * (1 + im*beta*tan(pi*alpha/2)*sign(k) * (abs(gamma*k)^(1-alpha) - 1)) + im*delta*k))
    else
        return(exp( -gamma*abs(k) * (1 + im*beta*(2/pi)*sign(k)*log(gamma*abs(k))) + im*delta*k))
}

int_ecf_dist_norm = function(x,mu,sigma)
{
    n=length(x)

    A=exp(-.25*x^2)%*%t(exp(-.25*x^2))
    lnB=.5*x%*%t(x)
    C=exp(log(A)+lnB)
    v = sum(C) * sqrt(pi)/(n*n) + sqrt(pi/(1+sigma*sigma)) - (2/n) * sqrt(pi/(1+sigma*sigma/2)) * sum(exp(-(x-mu)^2 / (4+2*sigma*sigma)))

    return(v)
}

reject_normality = function(x,n_rep=1000)
{

    ecf_x = ecf_fn(x)

    tuning = 10

    fweight = function(k){exp(-tuning*abs(k))}

    stable_max = optim(c(2*rbeta(1,1,1),rgamma(1,1,1)),function(pars) {int_ecf_dist(ecf=ecf_x,tcf=cstable,wtf=weight,alpha=pars[1],beta=0,delta=0,gamma=pars[2])$value},gr=NULL, method="L-BFGS-B", lower = c(.00001,.00001), upper = c(1.9999, 1000),control=list(factr=1))

    norm_max = optim(c(rgamma(1,1,1)),function(pars) {int_ecf_dist(ecf=ecf_x,tcf=cstable,wtf=weight,alpha=2.0,beta=0,delta=0,gamma=pars[1])$value},gr=NULL, method="L-BFGS-B", lower = c(.00001), upper = c(1000))

    sim_norm = matrix(rnorm(length(x)*n_rep,mean=0.0,norm_max$par[1]*sqrt(2)),nrow=n_rep,ncol=length(x))

    fit_norm = apply(sim_norm,1,function(y) {int_ecf_dist(ecf_fn(y),tcf=cstable,wtf=weight,alpha=2.0,beta=0,delta=0,gamma=norm_max$par[1])$value})


    null_cdf = ecdf(fit_norm)
    p_value = 1-null_cdf(norm_max$value)
    print(c(stable_max$value,norm_max$value,p_value,stable_max$par))
    return(list(pars=stable_max$par,p=p_value))

}

fpkm_test = function(fn)
{
    x = read.table(fn,header=F,row.names=1)
    x = log(x)
    #x = x[-which(rowSums(x)==-Inf),]
    x = x[apply(x,1,function(x){all(x>=1)}),]

    x = x - apply(x,1,median)
    res = apply(x[1:50,],1,reject_normality)
    return(res)
}

fpkm_test_dstable = function(fn)
{
    x = read.table(fn,header=F,row.names=1)
    x = log(x)
    x = x[apply(x,1,function(x){all(x>=1)}),]

    x = x - apply(x,1,median)
    resStable = apply(x[1:50,],1,function(dat) { optim(c(1,1),function(pars){-sum(dstable(dat,beta=0,alpha=pars[1],gamma=pars[2],delta=0,log=T))},gr=NULL,method="L-BFGS-B",lower=1e-5,upper=c(1.99999,1000))})

}

####################################
# fitting alpha-stable model (old) #
####################################

batch_optim = function(x,n_optim=10,center_x=TRUE)
{
    results = list()

    if (!is.matrix(x)) {
        x = t(as.matrix(x))
    }

    for (i in 1:nrow(x))
    {
        results[[i]] = list()
        results[[i]]$normal = list()
        results[[i]]$stable = list()
        if (center_x)
            x[i,] = x[i,] - median(x[i,])

        best_normal_idx = 1
        best_stable_idx = 1
        for (j in 1:n_optim)
        {
            print(c(i,j,best_normal_idx,best_stable_idx))
            optim_normal = optim(c(rexp(1)),function(pars){-sum(dstable(x[i,],beta=0,alpha=2,gamma=pars[1],delta=0,log=T))},gr=NULL,method="L-BFGS-B",lower=1e-5,upper=1000)

            optim_stable = optim(c(runif(1)*2, rexp(1)),function(pars){-sum(dstable(x[i,],beta=0,alpha=pars[1],gamma=pars[2],delta=0,log=T))},gr=NULL,method="L-BFGS-B",lower=c(1e-5,1e-5),upper=c(1.99999,1000))

            results[[i]]$normal[[j]] = optim_normal
            results[[i]]$stable[[j]] = optim_stable

            if (j > 1)
            {
                if (results[[i]]$normal[[best_normal_idx]]$value > optim_normal$value)
                    best_normal_idx = j

                if (results[[i]]$stable[[best_stable_idx]]$value > optim_stable$value)
                    best_stable_idx = j
            }
        }
        results[[i]]$best_normal_idx = best_normal_idx
        results[[i]]$best_stable_idx = best_stable_idx
        print(c(results[[i]]$best_normal_idx,results[[i]]$best_stable_idx))
    }

    return(results)
}

optim_fun_stable = function(dat) {
    optim_stable = optim(c(runif(1)*2, rexp(1),rnorm(1)),function(pars){-sum(dstable(dat,beta=0,alpha=pars[1],gamma=pars[2],delta=pars[3],log=T))},gr=NULL,method="L-BFGS-B",lower=c(1e-5,1e-5,-1000),upper=c(1.99999,1000,1000))
    return(optim_stable)
}

optim_fun_normal = function(dat) {
    optim_normal = optim(c(rexp(1),rnorm(1)),function(pars){-sum(dstable(dat,beta=0,alpha=2,gamma=pars[1],delta=pars[2],log=T))},gr=NULL,method="L-BFGS-B",lower=c(1e-5,-1000),upper=c(1000,1000))
    return(optim_normal)
}

batch_optim_mc = function(x,n_optim=10,num_core=1,center_x=TRUE)
{
    results = list()

    if (!is.matrix(x)) {
        x = t(as.matrix(x))
    }

    for (i in 1:nrow(x))
    {
        results[[i]] = list()

        if (center_x)
            x[i,] = x[i,] - median(x[i,])
        results[[i]]$data = x[i,]
        results[[i]]$name = rownames(x)[i]
        best_normal_idx = 1
        best_stable_idx = 1
        results[[i]]$normal = mclapply(1:n_optim,function(attempt) {optim_fun_normal(x[i,])},mc.cores=num_core,mc.preschedule=FALSE,mc.cleanup=9)
        results[[i]]$stable = mclapply(1:n_optim,function(attempt){optim_fun_stable(x[i,])},mc.cores=num_core,mc.preschedule=FALSE,mc.cleanup=9)
        for (j in 2:n_optim) {
            if (results[[i]]$normal[[best_normal_idx]]$value > results[[i]]$normal[[j]]$value && results[[i]]$normal[[j]]$convergence==0) {
                best_normal_idx = j
            }
            if (results[[i]]$stable[[best_stable_idx]]$value > results[[i]]$stable[[j]]$value && results[[i]]$stable[[j]]$convergence==0) {
                best_stable_idx = j
            }
        }

        results[[i]]$best_normal_idx = best_normal_idx
        results[[i]]$best_stable_idx = best_stable_idx
        print(c(i,results[[i]]$best_normal_idx,results[[i]]$best_stable_idx))
    }

    return(results)
}

getLikeRatio = function(curOptim) {
    bestStable = curOptim$best_stable_idx
    bestNormal = curOptim$best_normal_idx
    llStable = curOptim$stable[[bestStable]]$value
    llNormal = curOptim$normal[[bestNormal]]$value
    return(2*(llNormal-llStable))
}


batch_bootstrap = function(x,optims,ncores, nboots = 1000) {

    if (nrow(x)!=length(optims)) {
        print("Different number of optim results and genes")
        return(0)
    }

    p_values = c()

    bootstrap = matrix(nrow=0,ncol=nboots)

    for (i in 1:nrow(x)) {
        print(paste("currently analyzing gene",rownames(x)[i],sep=" "))
        curBestStable = optims[[i]]$best_stable_idx
        curBestNormal = optims[[i]]$best_normal_idx
        curSD = sqrt(2)*optims[[i]]$normal[[curBestNormal]]$par[1]
        curMu = optims[[i]]$normal[[curBestNormal]]$par[2]
        curAlpha = optims[[i]]$stable[[curBestStable]]$par[1]
        if (curAlpha >= 1.99) {
            p_values = c(p_values,1)
            bootstrap = rbind(bootstrap,rep(0,nboots))
            next
        }
        curFakeDat = matrix(nrow=nboots,ncol=ncol(x),byrow=T,rnorm(ncol(x)*nboots,mean=curMu,sd=curSD))
        curOptims = mclapply(1:nboots,function(curBoot) {batch_optim_mc(curFakeDat[curBoot,],n_optim=2,num_core=1,center_x=FALSE)},mc.cores=ncores,mc.preschedule=FALSE,mc.cleanup=9)
        curOptims = lapply(curOptims,function(x){x[[1]]})
        bootstrapLikeRatios = unlist(lapply(curOptims, getLikeRatio))
        curLikeRatio = getLikeRatio(optims[[i]])
        p_value = 1-ecdf(bootstrapLikeRatios)(curLikeRatio)
        p_values = c(p_values,p_value)
        bootstrap = rbind(bootstrap,bootstrapLikeRatios)
    }
    return(list(p=p_values,boot=bootstrap,dat=x))
}

sim_qts_old = function(nrep=50,nloc=8,nsam=8,popsize=2000,theta=.1,kernel=rnorm,scale=1,alpha=2.0)
{
    # loci
    b = list()
    for (n in 1:nloc)
    {
        b[[n]] = list()
        # samples
        for (m in 1:(nsam-1))
        {
            scale_per_locus = scale/((2^n)^(1/alpha))
            if (identical(kernel,rnorm))
                b[[n]][[m]] = batch_ms(nrep=nrep,s,pop=c(0,2^(m+1)),nloci=2^n,theta=theta,kernel=kernel,sd=scale_per_locus)

            else if (identical(kernel,rstable))
                b[[n]][[m]] = batch_ms(nrep=nrep,spop=c(0,2^(m+1)),nloci=2^n,theta=theta,kernel=kernel,gamma=scale_per_locus,alpha=alpha,beta=0)
            # get some traits from b
            # replicates
            for (l in 1:nrep)
            {
                cur_traits = b[[n]][[m]][[l]]$trait
                cur_traits = cur_traits-cur_traits[1,]
                cur_traits = cur_traits[2:2^(m+1),]
                print(cur_traits)
                # get KS from traits
                if (identical(kernel,rnorm))
                    ks_res = ks.test(cur_traits-median(cur_traits),pnorm,sd=sqrt(1/2*theta*scale^2))
                else if (identical(kernel,rstable))
                    ks_res = ks.test(cur_traits-median(cur_traits),pstable,gamma=scale*theta/2,alpha=alpha,beta=0)
                print(ks_res)
                b[[n]][[m]][[l]]$ks_stat = ks_res$statistic
                b[[n]][[m]][[l]]$ks_pvalue = ks_res$p.value
            }
        }
    }
    ks_stat=array(0,dim=c(nloc,nsam-1,nrep),dimnames=c("nloc","nsam","nrep"))
    ks_pvalue=array(0,dim=c(nloc,nsam-1,nrep),dimnames=c("nloc","nsam","nrep"))
    for (n in 1:nloc)
    {
        for (m in 1:(nsam-1))
        {
            for (l in 1:nrep)
            {
                ks_stat[n,m,l] = b[[n]][[m]][[l]]$ks_stat
                ks_pvalue[n,m,l] = b[[n]][[m]][[l]]$ks_pvalue
            }
        }
    }
    b$ks_stat = ks_stat
    b$ks_pvalue = ks_pvalue
    b$ks_stat_mean = apply(ks_stat,1:2,mean)
    b$ks_pvalue_mean = apply(ks_pvalue,1:2,mean)
    return(b)
}

#n is a vector of sample sizes at times corresponding to times in vector t
serial_coal = function (n, t, br = "coalescent", ...)
{
    total_lineages <- sum(n)
    nbr <- 2 * total_lineages - 2
    edge <- matrix(NA, nbr, 2)
    cur_sample_index = 1
    cur_time = t[cur_sample_index]
    cur_num_lineages = sum(n[1:cur_sample_index])
    pool = 1:cur_num_lineages
    num_lineages = sum(n[cur_sample_index:length(n)])
    nextnode = nbr + 1
    if (length(t) > 1) {
        next_time = t[cur_sample_index+1]
    } else {
        next_time = Inf
    }
    edge.length = numeric(nbr)
    i = 1
    h = numeric(nbr)
    h[1:n[1]] = t[1]
    while (cur_num_lineages > 1 | cur_time < t[length(t)]) {
        if (cur_num_lineages > 1) {
            #there is a possible coalescence
            coal_time = cur_time + rexp(1,choose(cur_num_lineages,2))
        } else {
            #otherwise, just pretend nothing happens
            coal_time = next_time + 1
        }
        if (coal_time < next_time) {
            #coalesce before new samples come in
            to_coal = sample(pool, size = 2)
            ind = (i-1)*2 + 1:2
            edge[ind,2] = to_coal
            edge[ind,1] = nextnode
            edge.length[ind] = coal_time - h[to_coal]
            cur_time = coal_time
            pool = c(pool[!pool%in%to_coal],nextnode)
            h[nextnode] = cur_time
            nextnode = nextnode - 1
            i = i + 1
            cur_num_lineages = cur_num_lineages - 1
            num_lineages = num_lineages - 1
        } else {
            #new samples come in
            new_samples = (sum(n[1:cur_sample_index])+1):(sum(n[1:cur_sample_index])+n[cur_sample_index+1])
            pool = c(pool,new_samples)
            cur_sample_index = cur_sample_index+1
            h[new_samples] = t[cur_sample_index]
            cur_time = t[cur_sample_index]
            if (length(t) > cur_sample_index) {
                next_time = t[cur_sample_index+1]
            } else {
                next_time = Inf
            }
            cur_num_lineages = cur_num_lineages + n[cur_sample_index]
        }
    }
    phy = list(edge = edge, edge.length = edge.length)
    phy$Nnode = total_lineages - 1
    phy$tip.label = 1:total_lineages
    class(phy) = "phylo"
    phy = reorder(phy)
    return(phy)
}

#This function requires the jumpdif.r for the function sim_data
sim_quantitative_trait = function(n,t,theta,num_loci,kernel,...) {
    #NB: theta = 4Nmu, mutation rate = 2Nmu!
    trait = numeric(sum(n))
    for (i in 1:num_loci) {
        cur_tree = serial_coal(n,t)
        cur_effect = sim_data(cur_tree,0,theta/2,kernel,...)
        trait = trait + cur_effect
    }
    return(trait)
}
