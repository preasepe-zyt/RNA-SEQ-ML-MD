#进行生存资料的特征筛选，从而进行特征基因提取
library(randomForestSRC)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(glmnet)
library(survival)
#字体
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)
train <- cbind(gsva_inflam$state.ch1, train[,-1])
names(train)[1] <- "group_list"
train$group_list <- as.factor(train$group_list)
#lasso回归
x1 <- data.matrix(train[, 2:ncol(train)])
x2 <- data.matrix(train[,1])
fit <- glmnet(x1, x2,
              family = "binomial",
              alpha = 1,
              relax=TRUE, 
              path=TRUE)
cvfit = cv.glmnet(x1, x2, 
                  nfold=10,
                  family = "binomial",
                  alpha=1,
                  type.measure = "deviance"
) 
#绘制lamba图
pdf("lambda_lasso.pdf",width = 10,height = 8,family="Arial")
plot(cvfit)
dev.off()

#用包自带的函数绘制lasso回归图
pdf("lasso.pdf",width = 10,height = 8,family="Arial")
plot(fit)
dev.off()
#最小值一个标准误时各变量的系数值
coef.min = coef(cvfit, s = "lambda.min") 
#非零系数
active.min = which(coef.min != 0)
#非零变量，基因名
geneids <- rownames(coef.min)[active.min]
#提取选中的基因对应的coefficient，回归系数
index.min = coef.min[active.min]
#按列合并两个变量
combine <- cbind(geneids, index.min)
#将数据转换为数据框
lasso.result<-as.data.frame(combine)
#将index.min数据类型改为数值型。
lasso.result$index.min <- as.numeric(lasso.result$index.min)
#将结果保存到文件中
write.table(lasso.result,"gene_index_lasso.txt",col.names=T,row.names=F, quote=F,sep="\t")

#ridge
x1 <- data.matrix(train[, 2:ncol(train)])
x2 <- data.matrix(train[,1])
fit2 <- glmnet(x1, x2,
              family = "binomial",
              relax=TRUE, 
              path=TRUE,
              alpha = 0)
cvfit2 = cv.glmnet(x1, x2, 
                  nfold=10,
                  family = "binomial",
                  alpha=0,
                  type.measure = "deviance"
) 
#绘制lamba图
pdf("lambda_ridge.pdf",width = 10,height = 8,family="Arial")
plot(cvfit2)
dev.off()

#用包自带的函数绘制ridge回归图
pdf("ridge.pdf",width = 10,height = 8,family="Arial")
plot(fit2)
dev.off()
#最小值一个标准误时各变量的系数值
coef.min2 = coef(cvfit2, s = "lambda.min") 
#非零系数
active.min2 = which(coef.min2 != 0)
#非零变量，基因名
geneids2 <- rownames(coef.min2)[active.min2]
#提取选中的基因对应的coefficient，回归系数
index.min2 = coef.min2[active.min2]
#按列合并两个变量
combine2 <- cbind(geneids2, index.min2)
#将数据转换为数据框
ridge.result <-as.data.frame(combine2)
#将index.min数据类型改为数值型。
ridge.result$index.min2 <- as.numeric(ridge.result$index.min2)
#将结果保存到文件中
write.table(ridge.result,"gene_index_ridge.txt",col.names=T,row.names=F, quote=F,sep="\t")

#elastic
fit3 <- glmnet(x1, x2,
               family = "binomial",
               relax=TRUE, 
               path=TRUE,
               alpha = 0.75)#根据0.5作为一个选择，看看偏向哪个
cvfit3 = cv.glmnet(x1, x2, 
                   nfold=10,
                   family = "binomial",
                   alpha=0.75,#根据0.5作为一个选择，看看偏向哪个
                   type.measure = "deviance"
) 
#绘制lamba图
pdf("lambda_elastic.pdf",width = 10,height = 8,family="Arial")
plot(cvfit3)
dev.off()


#用包自带的函数绘制elastic回归图
pdf("elastic.pdf",width = 10,height = 8,family="Arial")
plot(fit3)
dev.off()

#最小值一个标准误时各变量的系数值
coef.min3 = coef(cvfit3, s = "lambda.min") 
#非零系数
active.min3 = which(coef.min3 != 0)
#非零变量，基因名
geneids3 <- rownames(coef.min3)[active.min3]
#提取选中的基因对应的coefficient，回归系数
index.min3 = coef.min3[active.min3]
#按列合并两个变量
combine3 <- cbind(geneids3, index.min3)
#将数据转换为数据框
elastic.result<-as.data.frame(combine3)
#将index.min数据类型改为数值型。
elastic.result$index.min3 <- as.numeric(elastic.result$index.min3)
#将结果保存到文件中
write.table(elastic.result,"gene_index_elastic.txt",col.names=T,row.names=F, quote=F,sep="\t")

#随机森林
# Helper packages
library(dplyr)    # for data wrangling
library(ggplot2)  # for awesome graphics
# Modeling packages
library(ranger)   # a c++ implementation of random forest 
library(h2o) 
# number of features
n_features <- length(setdiff(names(train), "group_list"))

# create hyperparameter grid
hyper_grid <- expand.grid(
    mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
    min.node.size = c(1, 3, 5, 10), 
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(.5, .63, .8),                       
    rmse = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
        formula         = group_list ~ ., 
        data            = train, 
        num.trees       = n_features * 10,
        mtry            = hyper_grid$mtry[i],
        min.node.size   = hyper_grid$min.node.size[i],
        replace         = hyper_grid$replace[i],
        sample.fraction = hyper_grid$sample.fraction[i],
        verbose         = FALSE,
        seed            = 123,
        respect.unordered.factors = 'order',
    )
    # export OOB error 
    hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

# assess top 10 models
hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)

fit <- rfsrc(group_list ~ ., data = train,
             ntree = 3000, 
             mtry = 2,#grid research
             nodesize = 5, #grid research
             splitrule = "auc", #classification
             importance=TRUE,
             proximity = T,
             forest = T,
             standardized="TRUE"
             )

pdf("forest.pdf", width = 8, height = 8,family="Arial")
plot(fit)
dev.off()
#筛选特征基因
rftop <- var.select(fit)
rftop2 <- data.frame(
    Feature=rftop$topvars,
    vimp=rftop$varselect[rftop$topvars,2],
    Depth=rftop$varselect[rftop$topvars,1])
write.table(rftop2,"gene_index_randomforest.txt",col.names=T,row.names=F, quote=F,sep="\t")

#GBM模型
library(gbm)  
library(ranger)  
# search grid
hyper_grid <- expand.grid(
    n.trees = c(50, 100, 150),
    interaction.depth = c(3, 5, 7),
    shrinkage = c(0.01, 0.1, 0.2)
)
"""
train <- train %>%
    mutate(gbm = case_when(train$group_list=="Control"~0,
                        train$group_list=="Depression"~1))
train <- train[,-10]
"""
# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
        formula         = group_list ~ ., 
        data            = train, 
        num.trees       = n_features * 10,
        #mtry            = hyper_grid$mtry[i],
        #min.node.size   = hyper_grid$min.node.size[i],
        #replace         = hyper_grid$replace[i],
        sample.fraction = 0.5,
        verbose         = FALSE,
        seed            = 123,
        respect.unordered.factors = 'order',
    )
    # export OOB error 
    hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

# assess top 10 models
hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)



# results
arrange(hyper_grid, rmse)
# run a basic GBM model

ames_gbm1 <- gbm(
    formula = xg_g  ~ .,
    data = train,
    distribution = "bernoulli",  # SSE loss function
    n.trees = 50,
    shrinkage = 0.1,
    interaction.depth = 3,
    n.minobsinnode = 1,
    cv.folds = 10
)

# find index for number trees with minimum CV error
best <- which.min(ames_gbm1$cv.error)

# get MSE and compute RMSE
sqrt(ames_gbm1$cv.error[best])

# plot error curve
gbm.perf(ames_gbm1, method = "cv")


#xgb模型
library(xgboost)
library(Ckmeans.1d.dp)
xg_g <- case_when(train[,1]=="Control"~0,train[,1]=="Depression"~1)
# 定义参数网格
# hyperparameter grid
hyper_grid <- expand.grid(
    eta = 0.01,
    max_depth = 3, 
    min_child_weight = 3,
    subsample = 0.5, 
    colsample_bytree = 0.5,
    gamma = c(0, 1, 10, 100, 1000),
    lambda = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
    alpha = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
    rmse = 0,          # a place to dump RMSE results
    trees = 0          # a place to dump required number of trees
)

# grid search
for(i in seq_len(nrow(hyper_grid))) {
    set.seed(123)
    m <- xgb.cv(
        data = as.matrix(train[,2:ncol(train)]),
        label = xg_g,
        nrounds = 4000,
        objective = "binary:logistic",
        early_stopping_rounds = 50, 
        nfold = 10,
        verbose = 0,
        params = list( 
            eta = hyper_grid$eta[i], 
            max_depth = hyper_grid$max_depth[i],
            min_child_weight = hyper_grid$min_child_weight[i],
            subsample = hyper_grid$subsample[i],
            colsample_bytree = hyper_grid$colsample_bytree[i],
            gamma = hyper_grid$gamma[i], 
            lambda = hyper_grid$lambda[i], 
            alpha = hyper_grid$alpha[i]
        ) 
    )
    hyper_grid$rmse[i] <- min(m$evaluation_log$test_rmse_mean)
    hyper_grid$trees[i] <- m$best_iteration
}

# results
hyper_grid %>%
    filter(rmse > 0) %>%
    arrange(rmse) %>%
    glimpse()

# optimal parameter list
params <- list(
    eta = 0.01,
    max_depth = 3,
    min_child_weight = 3,
    subsample = 0.5,
    colsample_bytree = 0.5
)

# train final model

bstSparse <- xgboost(data = as.matrix(train[,2:ncol(train)]),
                     label = xg_g,
                     nrounds = 3944,
                     objective = "binary:logistic",
                     params = params,
                     verbose = 0)

#Perform the prediction
#重要重要性排序 
importance_matrix <- xgb.importance(model = bstSparse)  
pdf("xgb.pdf", width = 8, height = 8,family="Arial")
xgb.ggplot.importance(importance_matrix)+theme_classic()
dev.off()
xlsx::write.xlsx(importance_matrix,"xboost.xlsx")
# 设置10折交叉验证
# 创建XGBoost模型
model <- xgboost(data = as.matrix(X), label = y, objective = "binary:logistic")  # 适用于分类问题，如果是回归问题，使用"reg:linear"

# 设置10折交叉验证
control <- trainControl(method = "cv", number = 10)
# 执行交叉验证
result <- train(x = as.matrix(train[,2:ncol(train)]), y = xg_g, method = "xgbTree", trControl = control)

# 绘制交叉验证的性能曲线
plot(result)


####superpc####
library(superpc)
data <- list(x= t(train[,2:ncol(train)]),y=train$group_list,featurenames=colnames(train[,2:ncol(train)]))
fit <- superpc.train(data = data,type =  "regression") #default
cv.fit <- superpc.cv(fit,
                     data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
fit.red <- superpc.predict.red(cv.fit, 
                               data, 
                               data,
                               .6)
superpc.listfeatures(data, 
                     cv.fit,  
                     fit.red, 
                     num.features=10)

###SVM###
library(e1071)
library(kernlab)
library(caret)
library(tidyverse)

#SVM-RFE
#两个超参sigma (Sigma) C (Cost)
Profile=rfe(x= train[,2:ncol(train)],
            y= as.numeric(as.factor(train$group_list)),
            sizes = c(seq(1,8,by=1)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#plot
pdf(file="SVM-RFE.pdf", width=6, height=5.5,family="Arial")
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标记数量
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()
#保留特征基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

#PLS logistic regression 
# Perform 10-fold CV on a PLS model tuning the number of PCs to 
# use as predictors
set.seed(123)
cv_model_pls <- train(
    group_list ~ ., 
    data = train, 
    method = "pls",
    family = "binomial",
    trControl = trainControl(method = "cv", number = 10),
    preProcess = c("zv", "center", "scale"),
    tuneLength = 16
)

# Model with lowest RMSE
cv_model_pls$bestTune

# results for model with lowest loss
cv_model_pls$results %>%
    dplyr::filter(ncomp == pull(cv_model_pls$bestTune))

# Plot cross-validated RMSE
ggplot(cv_model_pls)
# Model interpretability packages
library(vip) 
re <- vip(cv_model_pls, num_features = 5)
write.table(file="pls-da.gene.txt", re$data, sep="\t", quote=F, row.names=F, col.names=F)
