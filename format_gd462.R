source("serial_coal.r")
require(ggplot2)
require(hexbin)
#require(GGally)

make_key = function(k_fn = "hgdp_key.txt")
{
    k = read.csv(k_fn,header=F,sep="\t",as.is=T)
    return(k)
}

#make_data = function(d_fn = "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt")
make_data = function(d_fn = "GD660.GeneQuantRPKM.txt", unique=T)
{
    d = read.csv(d_fn, header=T, sep="\t",as.is=T)
    r = d[,1]
    d = d[,-c(1,2,3,4)]
    n = colnames(d)
    if (unique)
    {
       n = unlist(lapply(n, function(x) { y=unlist(strsplit(x,'.',fixed=T))[1];return(y) }))
    }
    #r = d$TargetID
    d = matrix(as.numeric(unlist(d[,1:ncol(d)])),ncol=ncol(d))
    colnames(d) = n
    rownames(d) = r
    #d = d[,!duplicated(n)]
    idx = which( sapply(n,function(x){grepl("\\.1\\.M",x)}) )
    print(idx)
    #colnames(d)[!idx]
    print(colnames(d)[!idx])
    d = d[,!idx]
    d = d[,1:ncol(d)]
    return(d)
}

make_dupe_data = function(d_fn = "GD660.GeneQuantRPKM.txt")
{
    # read data
    d = read.csv(d_fn, header=T, sep="\t",as.is=T)
    r = d[,1]
    d = d[,-c(1,2,3,4)]
    n = colnames(d)

    # format individual ID
    n = unlist(lapply(n, function(x) { y=unlist(strsplit(x,'.',fixed=T))[1];return(y) }))

    # take difference of pairs of two
    dupecount = match(n, n)
    dupetable = as.data.frame( table(dupecount) )
    dupeidx = as.numeric(as.character( dupetable[dupetable$Freq > 1,]$dupecount ) )
    dupenames = n[dupeidx]
    y = matrix(0,ncol=length(dupenames),nrow=nrow(d))
    for (i in 1:10) { # length(dupeidx)) {

        j = as.integer(dupeidx[i])
        y[1:10,i] = log(d[1:10,j]) - log(d[1:10,j+1])
        print(c(i,j,dupeidx[i],colnames(d)[j], colnames(d)[j+1]))
        print(y[1:10,i])
        print(d[1:10,j])
        print(d[1:10,j+1])
        #u = rbinom(nrow(y),1,.5)*2 - 1
        #y[,i] = u * y[,i]
    }
    colnames(y) = dupenames
    rownames(y) = r

    ret = list()
    ret$y = y
    ret$d = d
    ret$dn = dupenames
    ret$di = dupeidx
    ret$dt = dupetable
    ret$dc = dupecount
    return(ret)

    return(y)
}

filter_data = function(d,cutoff=-5,f=1.0,filter=TRUE)
{
    d[,2:ncol(d)] = log(d[,2:ncol(d)])
    n = ncol(d)-1
    if (filter)
        above_cutoff = !apply( d[,2:ncol(d)],1,function(x) { sum(x>cutoff)/n < f })
    else {
        above_cutoff = !apply( d[,2:ncol(d)],1,function(x) { sum(x > -Inf)/n < f })
        #1:nrow(d)
        d[d == -Inf] = NA
    }
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
    d_Qn = apply(d,1,function(x) { Qn(x[!is.na(x)]) } )
    d_fn = apply(d,1,function(x) { f(x[!is.na(x)]) } )

    # get traits per Qn quantile
    qnt_val = quantile(d_Qn,probs=seq(0.0,1.0,0.1))
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
    #boxplot(new_d)
    plot(d_Qn, d_fn,col="black",xlab=expression(Q[n]),ylab=expression(KR[2]),cex=.5)
    for (i in 2:(n_bins-1)) {
        abline(v=qnt_val[i],lty=3)
    }
    abline(h=0.0,lty=2)
    return(ret)
}

qn_quant_plot = function(d, f=sk2, ylab, xlab=expression(Q[n])) {

    ret = list()

    # moments from data
    d_Qn = apply(d,1,function(x) { Qn(x[!is.na(x)]) } )
    d_fn = apply(d,1,function(x) { f(x[!is.na(x)]) } )

    # get traits per Qn quantile
    qnt_val = quantile(d_Qn,probs=seq(0.0,1.0,0.1))
    n_bins = length(qnt_val)
    qnt_elt = list()
    for (i in 1:(n_bins-1)) {
       x = intersect(which(d_Qn > qnt_val[i]), which(d_Qn <= qnt_val[i+1]))
       qnt_elt[[i]] = x
    }

    dev.new()
    df = data.frame(x=d_Qn,y=d_fn)
    p = ggplot(df, aes(x=x,y=y)) + stat_binhex(bins=100) + geom_vline(xintercept=qnt_val[2:(n_bins-1)], col="black", linetype="dotted") + geom_hline(yintercept=median(d_fn),col="black",linetype="dashed") + xlab(xlab) + ylab(ylab) + scale_fill_gradient(low="gray",high="black")
    print(p)
    return(ret)
}

quant_kr_sk = function(d,th=.5,dq=0.02) {

    ret = list()

    # moments from data
    d_Qn = apply(d,1,function(x) { Qn(x[!is.na(x)]) } )
    d_sk = abs(apply(d,1,function(x) { sk2(x[!is.na(x)]) } ))
    d_kr = abs(apply(d,1,function(x) { kr2(x[!is.na(x)]) } ))

    # get traits per Qn quantile
    qnt_qn = quantile(d_Qn,probs=seq(0.0,1.0,dq))
    qnt_sk = quantile(d_sk,probs=th)
    qnt_kr = quantile(d_kr,probs=th)

    # get skew/kurt upper 50 per bracket
    n_bins = length(qnt_qn)
    qnt_elt = matrix(c(0),n_bins-1,4)
    for (i in 1:(n_bins-1)) {
        x = intersect(which(d_Qn > qnt_qn[i]), which(d_Qn <= qnt_qn[i+1]))
        tmp_sk = d_sk[x]
        tmp_kr = d_kr[x]
        lbl = sapply(1:length(x), function(y) {
        if (tmp_sk[y] < qnt_sk && tmp_kr[y] < qnt_kr) {
            return(1)
        } else if (tmp_sk[y] < qnt_sk && tmp_kr[y] >= qnt_kr) {
            return(2)
        } else if (tmp_sk[y] >= qnt_sk && tmp_kr[y] < qnt_kr) {
            return(3)
        } else if (tmp_sk[y] >= qnt_sk && tmp_kr[y] >= qnt_kr) {
            return(4)
        }})
        v = sapply(1:4, function(y) { sum(lbl==y) })
        qnt_elt[i,] = v
    }
    v = data.frame(quantile=(1:(n_bins-1))*dq*100-1,s0k0=qnt_elt[,1], s0k1=qnt_elt[,2], s1k0=qnt_elt[,3], s1k1=qnt_elt[,4])
    mv = melt(v,id=c("quantile"))
    
    #txt_lbls = c(   paste("S<",th,",K<",th,sep=""),paste("S<",th,",K>",th,sep=""),paste("S>",th,",K<",th,sep=""),paste("S>",th,",K>",th,sep=""))
    #txt_lbls = c("S-K-","S-K+","S+K-","S+K+")
    txt_lbls = c("symmetric, mesokurtic","symmetric, non-mesokurtic","skewed, mesokurtic","skewed, non-mesokurtic")

    mn_qnt_qn_txt = as.character(sapply(2:n_bins, function(y) { mean(qnt_qn[(y-1):y]) } ))
    p = ggplot(mv,aes(x=quantile,y=value,fill=variable)) +
        geom_area(position="stack") +
        scale_fill_discrete(name="",labels=txt_lbls) +
        scale_x_continuous(breaks=seq(0,100,10)) +
        xlab(expression(Q[n] ~ quantiles)) +
        ylab(No. ~ SK[2] ~ and ~ KR[2] ~ matches) +
        ggtitle(paste(th*100,"th skewness and kurtosis quantiles",sep="")) +
        geom_hline(yintercept=c(th*th*sum(v[1,])),linetype="dashed") +
        guides(fill=guide_legend(reverse=T))
    print(p)

    return(v)
}
# d: data frame, traits x individuals
# f: function you want to compare vs Qn quantiles -- e.g. sk2, kr2
# g: function you want to use to summarize f-quantile results -- e.g. var, median
q_var_moment2 = function(d,dq=0.02) {

    # moments from data
    d_Qn = apply(d,1,function(x) { Qn(x[!is.na(x)]) } )
    d_sk = apply(d,1,function(x) { sk2(x[!is.na(x)]) } )
    d_kr = apply(d,1,function(x) { kr2(x[!is.na(x)]) } )

    # sim data
    dr = 1
    sim = data.frame( matrix( 0, nrow(d)/dr, ncol(d) ) )
    for (i in 1:(nrow(d)/dr)) {
        sim[i,] = rnorm(ncol(d), 0., d_Qn[i])
    }
    d_sk_null = apply(sim,1,function(x) { sk2(x) } )
    d_kr_null = apply(sim,1,function(x) { kr2(x) } )

    # get traits per Qn quantile
    qnt_val = quantile(d_Qn,probs=seq(0.0,1.0,dq))
    n_bins = length(qnt_val)
    qnt_elt = matrix(c(0),n_bins-1,5)
    for (i in 1:(n_bins-1)) {
        x = intersect(which(d_Qn > qnt_val[i]), which(d_Qn <= qnt_val[i+1]))
        yy = as.integer(x/dr)
        qnt_elt[i,] = c(i*dq*100-1,var(d_sk[x]), var(d_kr[x]), var(d_sk_null[yy]), var(d_kr_null[yy]))
    }

    df = data.frame(    n=qnt_elt[,1],
                        sk=qnt_elt[,2],
                        kr=qnt_elt[,3],
                        skn=qnt_elt[,4],
                        krn=qnt_elt[,5])

    p = ggplot(df, aes(x=n)) + 
        geom_line(aes(y=sk,colour="sk2",linetype="eur")) + 
        geom_line(aes(y=kr,colour="kr2",linetype="eur")) + 
        geom_line(aes(y=skn,colour="sk2",linetype="sim")) +
        geom_line(aes(y=krn,colour="kr2",linetype="sim")) +
        scale_colour_discrete(name="Moment", breaks=c("sk2","kr2"), labels=c(expression(SK[2]),expression(KR[2]))) + 
        scale_linetype_discrete(name="Data",breaks=c("eur","sim"),labels=c("empirical","simulated")) + 
        scale_x_continuous(breaks=seq(0,100,10)) +
        xlab(expression(Q[n] ~ quantiles)) +
        ylab("Sample moment variance")
    print(p)

    return(qnt_elt)
}
#
## workflow
if (FALSE) 
{
# make data
k=make_key()
d=make_data()
df=filter_data(d)
eur=make_pop_data(k,df,not=T)
eur=eur[,-1] # crazy outlier
yri=make_pop_data(k,df)

# SW discoveries
eur_sw_fdr_idx = get_sw_sig(eur, qval=0.05, ndrop=0, ns=ncol(eur))
yri_sw_fdr_idx = get_sw_sig(yri, qval=0.05, ndrop=0)
eurss_sw_fdr_idx = get_sw_sig(eur, qval=0.05, ndrop=, ns=ncol(yri))

# naive robust moments
eur_sk = apply(eur,1,sk2)
eur_kr = apply(eur,1,kr2)
eur_qn = apply(eur,1,Qn)
#plot(eur_qn,eur_sk,col="gray")

# Wilcox test
#hist(eur_sk,breaks=100,xlab=expression("SK"[2]),main="");abline(v=0,lw=3)
#hist(eur_kr,breaks=100,xlab=expression("KR"[2]),main="",xlim=c(-1,4));abline(v=0,lw=3)

# scatter-quantile w/ hexbins
sq_sk2 = qn_quant_plot(eur, f=sk2, ylab=expression(SK[2]))
sq_kr2 = qn_quant_plot(eur, f=kr2, ylab=expression(KR[2]))

# stacked histograms of gene expression quantiles
qks_0_5 = quant_kr_sk(eur, th=0.5)
qks_0_95 = quant_kr_sk(eur, th=0.95)

# quantile for var of moments
qvm2 = q_var_moment2(eur)

# check w/in variation
yy = make_dupe_data()[rownames(df),]
fy = fy[,-1] # drop outlier, 168 pairs
fy = apply(fy,2,function(x) { x/eur_qn }) # divide each within-ind trait by among-ind variance
kry = apply(fy,2,kr2) # kurtosis of within-ind traits for each individual
sky = apply(fy,2,sk2) # skewness of within-ind traits for each individual


}


## quantile summary stats
#eur_sk_med = q_var_moment(eur, f=sk2, g=median)
#eur_sk_var = q_var_moment(eur, f=sk2, g=var)
#eur_kr_med = q_var_moment(eur, f=kr2, g=median)
#eur_kr_var = q_var_moment(eur, f=kr2, g=var)
#
## index ordered estimates
#order_sk_idx = order(eur_sk)
#order_kr_idx = order(eur_kr)
#
#plot(eur_kr_var$d, col=1, type="l", ylim=c(-0.1,0.1))
#lines(eur_kr_var$sim, col=1, lty=2)
#lines(eur_sk_var$d, col=2, lty=2)
#lines(eur_sk_var$sim, col=2, lty=2)
#lines(eur_kr_med$d, col=3, lty=2)
#lines(eur_kr_med$sim, col=3, lty=2)
#lines(eur_sk_med$d, col=4, lty=2)
#lines(eur_sk_med$sim, col=4, lty=2)
#
## dip discoveries
#eur_dip=get_dip_pvals(eur,drop_outlier=F,ncores=4)
#eur_dip_fdr = fdrtool(0-eur_dip$dip, statistic="pvalue", cutoff.method="fndr")
#eur_dip2 = c(dipp.tantrum(eur,dip(eur),M=nrow(eur)*1.5)$p.value)
#
#
## compare kurt
#
#
#
## compare GWAS
#
#
## robust
#
#
## TODO
## 1) SW test, rampant non-normality
## 2) skew/kurt wrt Qn, increases
## 3) multimodality, dip test
## 4) re-check GWAS+dip enrichment
## 5) example histograms
#
#
#
#}
#
#
#
##eur_sw = list()
##for (i in 6:6) { eur_sw[[i]] = get_sw_pvals(eur,T,n=i) }
##lapply(eur_sw, function(x) { sum(x<0.005) })
##yri_sw = list()
##for (i in 10:10) { yri_sw[[i]] = get_sw_pvals(yri,T,n=i) }
##lapply(yri_sw, function(x) { sum(x<0.005) })
## Drop most outliers as possible
## idx 7 has 4700 discoveries
## idx 10 has 636 discoveries
#
##eur_sw_fdr = fdrtool(eur_sw[[6]]$sw,statistic=c("pvalue"),cutoff.method=c("fndr"))
##yri_sw_fdr = fdrtool(yri_sw[[10]]$sw,statistic=c("pvalue"),cutoff.method=c("fndr"))
##eur_sw_fdr_idx = which( eur_sw_fdr$qval < 0.005)
##yri_sw_fdr_idx = which( yri_sw_fdr$qval < 0.005)
