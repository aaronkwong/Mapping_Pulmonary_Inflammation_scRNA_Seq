##############################################
#use R 4.0.3
##############################################

library(renv)

options(install.opts = "--no-lock")

renv::init(bare=TRUE)
# "y" answer yes to continue prompts
renv::restore(exclude=c("nebula","Augur","CellChat","DoubletFinder","sctransform","BPSC","MatrixGenerics","presto","scClustViz","sparseMatrixStats"),prompt=FALSE)
# y
# y

library(devtools)
install_github("lhe17/nebula@cb8fd6622825d325acf89694102f1f44d1186b76")
 
install_github("chris-mcginnis-ucsf/DoubletFinder@5dfd96b06365d7843adf3a72ffb6a30f42c74a01")
 
install_github("sqjin/CellChat@af7e15246f29904a99047d74c1f6126d0f831bbc")
 
install_github("neurorestore/Augur@8539dd1ddaf44e245b97a88256a77093e816fbe7")

install_github("BaderLab/scClustViz@32064ce623fb0388e7647fe0ede52d0dab3036ab")