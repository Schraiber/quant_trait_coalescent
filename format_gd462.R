require(ggplot2)

make_key = function(k_fn = "hgdp_key.txt")
{
    k = read.csv(k_fn,header=F,sep="\t",as.is=T)
    return(k)
}

make_data = function(d_fn = "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt")
{
    d = read.csv(d_fn, header=T, sep="\t",as.is=T)
    n = colnames(d)[5:ncol(d)]
    r = d$TargetID
    d = matrix(as.numeric(unlist(d[,5:ncol(d)])),ncol=ncol(d)-4)
    colnames(d) = n
    rownames(d) = r
    return(d)
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

hist_compare = function(d1,d2,names=c("YRI","EUR"))
{
    df1=data.frame(gexp=d1)
    df2=data.frame(gexp=d2)
    df1$pop=names[1]
    df2$pop=names[2]
    gexp=rbind(df1,df2)
    #ggplot(gexp,aes(gexp,fill=pop))+geom_bar(pos="dodge")
    ggplot(gexp,aes(gexp,fill=pop))+geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
}

# just some basic ideas
hmm = function(d1,d2)
{
    m1_2=get_expr_var(d1)
    m1_3=get_expr_skew(d1)
    m1_4=get_expr_kurt(d1)
    dip1=get_dip_pvals(d1,F,F)
    sw1=get_sw_pvals(d1)
   
    m2_2=get_expr_var(d2)
    m2_3=get_expr_skew(d2)
    m2_4=get_expr_kurt(d2)
    dip2=get_dip_pvals(d2,F,F)
    sw2=get_sw_pvals(d2)

    dev.new();plot(m1_2[,1],m2_2[,1],xlab="EUR",ylab="YRI",main="var")
    dev.new();plot(m1_3[,1],m2_3[,1],xlab="EUR",ylab="YRI",main="skew")
    dev.new();plot(m1_4[,1],m2_4[,1],xlab="EUR",ylab="YRI",main="kurt")

    dev.new();plot(dip1[,1],dip2[,1],xlab="EUR",ylab="YRI",main="dip")
    dev.new();plot(sw1[,1],sw2[,1],xlab="EUR",ylab="YRI",main="SW")

}
