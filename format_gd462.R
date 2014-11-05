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
