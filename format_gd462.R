source("serial_coal.r")
require(ggplot2)
require(GGally)

make_key = function(k_fn = "hgdp_key.txt")
{
    k = read.csv(k_fn,header=F,sep="\t",as.is=T)
    return(k)
}

#make_data = function(d_fn = "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt")
make_data = function(d_fn = "GD660.GeneQuantRPKM.txt", unique=T)
{
    d = read.csv(d_fn, header=T, sep="\t",as.is=T)
    d = d[,-c(2,3,4)]
    n = colnames(d)
    if (unique)
    {
       n = unlist(lapply(n, function(x) { y=unlist(strsplit(x,'.',fixed=T))[1];return(y) }))
    }
    r = d$TargetID
    d = matrix(as.numeric(unlist(d[,1:ncol(d)])),ncol=ncol(d))
    colnames(d) = n
    rownames(d) = r
    d = d[,!duplicated(n)]
    d = d[,2:ncol(d)]
    return(d)
}

filter_data = function(d,cutoff=-5,f=1.0)
{
    d[,2:ncol(d)] = log(d[,2:ncol(d)])
    n = ncol(d)-1
    above_cutoff = !apply( d[,2:ncol(d)],1,function(x) { sum(x>cutoff)/n < f })
    return(d[above_cutoff,])
}

# CEU, FIN, GBR, TSI, YRI
match_cols = function(k,d)
{
    n = colnames(d)
    m = n
    for (i in 1:length(n))
    {
        idx = unlist(k[which(k[,2]==n[i]),1])
        if (is.na(m[i]) || identical(idx,character(0)))
            m[i] = n[i]
        else
            m[i] = idx
    }
    return(m)
}

make_pop_data = function(k,d,s="YRI",not=F)
{
    idx = c()
    idx = which(match_cols(k,d)==s)
    if (not)
        idx = -idx
    return(d[,idx])
}

hist_compare = function(d1,d2,names=c("EUR","YRI"))
{
    df1=data.frame(gexp=d1)
    df2=data.frame(gexp=d2)
    df1$pop=names[1]
    df2$pop=names[2]
    gexp=rbind(df1,df2)
    #ggplot(gexp,aes(gexp,fill=pop))+geom_bar(pos="dodge")
    #ggplot(gexp,aes(gexp,fill=pop))+geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
    ggplot(gexp,aes(gexp,fill=pop))+geom_histogram(alpha = 0.5, position = 'identity')
}

get_sw_sig = function(d, qval=0.005, ndrop=0, ns=0) {
    if (ns== 0) ns = ncol(d)
    tmp = matrix(nrow=nrow(d),ncol=ns)
    for (i in 1:nrow(d)) {
       tmp[i,] = sample(d[i,],ns) 
    }
    sw = get_sw_pvals(tmp, T, ndrop)
    fdr = fdrtool(sw$sw, statistic="pvalue", cutoff.method="fndr")
    idx = which(fdr$qval < qval)
}

# just some basic ideas
hmm = function(d1,d2)
{
    m1_1=apply(d1[,1:ncol(d1)],1,function(x){mean(x,na.rm=T)})
    m2_1=apply(d2[,1:ncol(d2)],1,function(x){mean(x,na.rm=T)})
    #dev.new();plot(m1_1,m2_1,xlab="EUR",ylab="YRI",main="mean")

    m1_2=get_expr_var(d1)
    m2_2=get_expr_var(d2)
    #dev.new();plot(m1_2[,1],m2_2[,1],xlab="EUR",ylab="YRI",main="var")

    m1_3=get_expr_skew(d1)
    m2_3=get_expr_skew(d2)
    #dev.new();plot(m1_3[,1],m2_3[,1],xlab="EUR",ylab="YRI",main="skew")

    m1_4=get_expr_kurt(d1)
    m2_4=get_expr_kurt(d2)
    #dev.new();plot(m1_4[,1],m2_4[,1],xlab="EUR",ylab="YRI",main="kurt")

    dip1=get_dip_pvals(d1,F,F)
    dip2=get_dip_pvals(d2,F,F)
    #dev.new();plot(1-dip1[,1],1-dip2[,1],xlab="EUR",ylab="YRI",main="dip")

    sw1=get_sw_pvals(d1)
    sw2=get_sw_pvals(d2)
    #dev.new();plot(sw1[,1],sw2[,1],xlab="EUR",ylab="YRI",main="SW")

    d = data.frame( 
                    mean1=m1_1,mean2=m2_1,
                    var1=m1_2[,1],var2=m2_2[,1],
                    skew1=m1_3[,1],skew2=m2_3[,1],
                    kurt1=m1_4[,1],kurt2=m2_4[,1],
                    dip1=1-dip1[,1],dip2=1-dip2[,1],
                    sw1=sw1[,1],sw2=sw2[,1]
                )
    dev.new()
    ggpairs(data=d,columns=1:ncol(d))
}

# d: data frame, traits x individuals
# f: function you want to compare vs Qn quantiles -- e.g. sk2, kr2
# g: function you want to use to summarize f-quantile results -- e.g. var, median
q_var_moment = function(d, f=sk2, g=var) {

    ret = list()

    # moments from data
    d_Qn = apply(d,1,Qn)
    d_fn = apply(d,1,f)

    # get traits per Qn quantile
    qnt_val = quantile(d_Qn,probs=seq(0.0,1.0,0.05))
    n_bins = length(qnt_val)
    qnt_elt = list()
    for (i in 1:(n_bins-1)) {
       x = intersect(which(d_Qn > qnt_val[i]), which(d_Qn <= qnt_val[i+1]))
       qnt_elt[[i]] = x
    }

    # null
    sim = c()
    for (i in 1:(n_bins-1)) {
        sigma = (qnt_val[i]+qnt_val[i+1])/2
        x = c()
        for (j in 1:length(qnt_elt[[i]])) { 
            x[j] = f(rnorm(ncol(d), 0., sigma))
        }
        sim[i] = g(x)
    }

    # get variance of
    ret$d = unlist(lapply(qnt_elt, function(x) { g(d_fn[x]) } ))
    ret$dU = unlist(lapply(qnt_elt, function(x) { quantile(d_fn[x], 0.7) }))
    ret$dL = unlist(lapply(qnt_elt, function(x) { quantile(d_fn[x], 0.3) }))
    ret$sim = sim

    new_d = lapply(qnt_elt, function(x) { d_fn[x] })
    dev.new()
    boxplot(new_d)
    return(ret)
}


# workflow
if (false)
{
# make data
k=make_key()
dd=make_data()
df=filter_data(dd)
eur=make_pop_data(k,df,not=T)
eur=eur[,-1] # crazy outlier
yri=make_pop_data(k,df)

# SW discoveries
eur_sw_fdr_idx = get_sw_sig(eur, qval=0.005, ndrop=6, ns=ncol(yri))
yri_sw_fdr_idx = get_sw_sig(yri, qval=0.005, ndrop=10)

# dip discoveries
eur_dip=get_dip_pvals(eur,drop_outlier=F,ncores=4)
eur_dip_fdr = fdrtool(0-eur_dip$dip, statistic="pvalue", cutoff.method="fndr")
eur_dip2 = c(dipp.tantrum(eur,dip(eur),M=nrow(eur)*1.5)$p.value)

# quantile stuff
eur_sk_med = q_var_moment(eur, f=sk2, g=median)
eur_sk_med = q_var_moment(eur, f=sk2, g=var)
eur_kr_var = q_var_moment(eur, f=kr2, g=median)
eur_kr_var = q_var_moment(eur, f=kr2, g=var)

plot(eur_kr_var$d, col=1, type="l", ylim=c(-0.1,0.1))
lines(eur_kr_var$sim, col=1, lty=2)
lines(eur_sk_var$d, col=2, lty=2)
lines(eur_sk_var$sim, col=2, lty=2)
lines(eur_kr_med$d, col=3, lty=2)
lines(eur_kr_med$sim, col=3, lty=2)
lines(eur_sk_med$d, col=4, lty=2)
lines(eur_sk_med$sim, col=4, lty=2)

# compare kurt



# compare GWAS


# robust
}



#eur_sw = list()
#for (i in 6:6) { eur_sw[[i]] = get_sw_pvals(eur,T,n=i) }
#lapply(eur_sw, function(x) { sum(x<0.005) })
#yri_sw = list()
#for (i in 10:10) { yri_sw[[i]] = get_sw_pvals(yri,T,n=i) }
#lapply(yri_sw, function(x) { sum(x<0.005) })
# Drop most outliers as possible
# idx 7 has 4700 discoveries
# idx 10 has 636 discoveries

#eur_sw_fdr = fdrtool(eur_sw[[6]]$sw,statistic=c("pvalue"),cutoff.method=c("fndr"))
#yri_sw_fdr = fdrtool(yri_sw[[10]]$sw,statistic=c("pvalue"),cutoff.method=c("fndr"))
#eur_sw_fdr_idx = which( eur_sw_fdr$qval < 0.005)
#yri_sw_fdr_idx = which( yri_sw_fdr$qval < 0.005)
