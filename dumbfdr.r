get_fdr = function(x,alpha=0.2) {
	x_sort=sort(x)
	fdr=NULL
	for (i in 1:length(x)) { 
		fdr=c(fdr,x[i] >= match(x[i],x_sort)*alpha/length(x))
	}
	return(fdr)
}


d1=read.table("p_vals.txt")
d2=read.table("p_vals_remove.txt")
p1=d1$V2
p2=d2$V2


get_fdr_mtx = function(p) {

	fdrs = matrix(0,nrow=21,ncol=2)
	for (i in 1:21) {
		alpha = (i-1)*0.05
		fdr1 = get_fdr(p,alpha)
		fdrs[i,1]=alpha
		fdrs[i,2]=sum(fdr1==FALSE)
	}
	colnames(fdrs)=c("alpha","n_discovery")
	
	return(fdrs)
}

get_gexp = function(fn="LA_normalizedExpression.txt",cutoff=-5.0)
{
	d = read.table(fn,row.names=1,header=FALSE)
	d[,1:ncol(d)] = log(d[,1:ncol(d)])
	above_cutoff = !apply(d[,1:ncol(d)], 1, function(x) { any(x <= cutoff) })
	return(d[above_cutoff,])
}

get_pvals = function(dd,drop_outlier=TRUE)
{
	if (drop_outlier)
		r = apply( dd[,1:ncol(dd)],1,function(x) { 
			shapiro.test(x[-which.max(abs(median(x)-x))])
		})
	else
		r = apply(dd[,1:ncol(dd)],1,shapiro.test)
	p = unlist(lapply(r,function(x){x$p.value}))
	df = data.frame(V1=dd[,1],V2=unlist(p))
	return(df)
}
