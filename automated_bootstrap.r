source("serial_coal.r")
library(compiler)
enableJIT(3)

args = commandArgs(trailingOnly=TRUE)

to_load = args[1]
output_destination = args[2]
num_cores = args[3]

load(to_load)
gene.boot = list()
if (length(gene.optim)) {
   n_genes = length(gene.optim)
   n_sam = length(gene.optim[[1]]$data)
   dat = unlist(lapply(gene.optim,function(x){x$data}))
   gene.dat = matrix(nrow=n_genes,ncol=n_sam,byrow=T,dat)
   rownames(gene.dat) = unlist(lapply(gene.optim,function(x){x$name}))
   #print(gene.dat)
   gene.boot = batch_bootstrap(gene.dat,gene.optim,num_cores)
}

save(gene.boot,file=output_destination)
quit()
