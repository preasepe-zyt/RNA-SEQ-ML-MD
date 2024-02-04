
# install packages 这三个安装不成功的话，就安后面的bseqsc包也行
install.packages('e1071')
install.pacakges('parallel')
install.packages('preprocessCore')
library(e1071)
library(preprocessCore)
library(parallel)
library(devtools)
library(bseqsc)#这个包携带大量CIBERSORT的依赖，前三个安装不好可以安装他

#启动
source("CIBERSORT.R") #启动这个函数，必须在哦那个一个文件夹内才可哟
results <- CIBERSORT(sig_matrix ="LM22.txt", mixture_file ="ids_exprs.txt", perm = 1000, QN = T)
