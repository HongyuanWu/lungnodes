
data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/lungnodes/master/extdata/TCGA.3108.Marker.DJG.beta.txt",head=T)
colnames(data)[2]<-"phen"
iid<-data[,1]
phen<-data[,2]
table(phen)
clinical<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/lungnodes/master/extdata/clinical.txt",head=T,row.names = 1)
input<-data.matrix(data[,2:ncol(data)])
input[input[,1]<3,1]<-0
input[input[,1]==3,1]<-1
input[1:5,1:5]
match(iid,rownames(clinical))

library("randomForest")
library("neuralnet")
library("arm")
library("PredictABEL")

set.seed(49)
cv.error <- NULL
k <- 2
rlt1<-c()
rlt2<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  train.cv[1:3,1:3]
  test.cv[1:3,1:3]
  
  P=apply(train.cv[,2:ncol(train.cv)],2,function(x) summary(glm(as.factor(train.cv[,1])~x,family=binomial))$coefficients[2,4])
  
  train.cv<-train.cv[,c(1,match(names(P[head(order(P),n=5200)]),colnames(train.cv)))]
  test.cv<-test.cv[,c(1,match(names(P[head(order(P),n=5200)]),colnames(test.cv)))]

  train.cv[1:5,1:5]
  
  RF <- randomForest(as.factor(phen) ~ ., data=data.frame(na.omit(train.cv)), importance=TRUE,proximity=T)
  imp<-RF$importance
  head(imp)
  imp<-imp[order(imp[,4],decreasing = T),]
  
  topvar<-match(rownames(imp)[1:30],colnames(input))
  
  train.cv <- input[index,c(1,topvar)]
  test.cv <- input[-index,c(1,topvar)]
  
  RF <- randomForest(as.factor(phen) ~ ., data=data.frame(na.omit(train.cv)), importance=TRUE,proximity=T)
  RF
  RF <- randomForest(as.factor(phen) ~ ., data=data.frame(na.omit(train.cv)), importance=TRUE,proximity=T)
  
  n <- colnames(train.cv)
  f <- as.formula(paste("phen ~", paste(n[!n %in% "phen"], collapse = " + ")))
  
  nn <- neuralnet(f,data=train.cv,hidden=c(5,3),act.fct = "logistic",linear.output = T)
  pr.nn <- neuralnet::compute(nn,test.cv)
  trainRlt<-data.frame(phen=train.cv[,1],pred=unlist(nn$net.result[[1]][,1]))
  testRlt<-data.frame(phen=test.cv[,1],pred=unlist(pr.nn$net.result[,1]))
  rownames(trainRlt)=row.names(train.cv)
  rownames(testRlt)=row.names(test.cv)
  rlt1<-rbind(rlt1,trainRlt)  
  rlt2<-rbind(rlt2,testRlt)
  print(i)
}

data1<-na.omit(data.frame(rlt1))
data2<-na.omit(data.frame(rlt2))
model.glm1 <- bayesglm(phen~.,data=rlt1,family=binomial(),na.action=na.omit)
model.glm2 <- bayesglm(phen~.,data=rlt2,family=binomial(),na.action=na.omit)
pred1 <- predRisk(model.glm1)
pred2 <- predRisk(model.glm2)
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))

### ROC 
pdf("mRNA.drugresponse.pdf")
par(mfrow=c(2,2),cex.lab=1.5,cex.axis=1.5)
plotROC(data=data1,cOutcome=1,predrisk=cbind(pred1))
plotROC(data=data2,cOutcome=1,predrisk=cbind(pred2))
dev.off()

### heatmap
library("gplots")
source("https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/HeatMap.R")
input<-na.omit(input)
temp<-t(input[,2:ncol(input)])
colnames(temp)<-input[,1]
pdf("heatmap.lungnodes.pdf")
HeatMap(temp)
dev.off()

chr<-unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"_"), function(x) x[1]))
pos<-unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"_"), function(x) x[2]))
probe<-data.frame(CHR=chr,POS=pos)
head(probe)
