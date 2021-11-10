# GSE130563 data cleaning

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

geo_id <- read.csv("C:\\Users\\Steven\\Desktop\\github\\LASSO\\GSE130563_gene id.csv", head = T)
geo_name <- bitr(geo_id[,], fromType = "ENTREZID", toType = "SYMBOL", OrgDb="org.Hs.eg.db", drop = F)
write.csv(geo_name, "C:\\Users\\Steven\\Desktop\\github\\LASSO\\GSE130563_gene name.csv")

gse130563 <- read.csv("C:\\Users\\Steven\\Desktop\\github\\LASSO\\GSE130563.csv", head = F)
gse130563 <- t(gse130563)
gene_data <- apply(gse130563[2:nrow(gse130563),3:ncol(gse130563)],1:2, as.double)

exp_2 <- function(x){
  2 ^ x
}

gene_data <- apply(gene_data,1:2, exp_2)
colnames(gene_data) <- t(geo_name[,2])
rownames(gene_data) <- gse130563[2:47,1]
gene_data <- as.data.frame(gene_data)


# volcano plot

library(ggplot2)
cachexia <- gene_data[17:33,]
noncachexia <- gene_data[c(1:16, 34:38),]

fold_change <- apply(cachexia, 2, mean) / apply(noncachexia, 2, mean)
log_fc <- log10(fold_change)

p_value <-c()
for(i in 1:ncol(gene_data)){
  p_value <- c(p_value, t.test(gene_data[17:33,i], gene_data[c(1:16, 34:38),i], var.equal = T)$p.value)
}
log_p <- (-1) * log10(p_value)
gene_data_2 <- rbind(log_fc, log_p)
gene_data_2 <- t(gene_data_2)

vol <- ggplot(as.data.frame(gene_data_2), aes(x = log_fc, y = log_p)) + 
  geom_point(color = ifelse(abs(log_fc) > log10(1.5), "blue", "black")) + 
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dotted") + 
  geom_vline(xintercept = c(log10(1.5), -log10(1.5)), col = "red", linetype = "dotted") + 
  labs(title = "Volcano plot") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), axis.title = element_text(size = 14)) + 
  xlab("log fold change") + 
  ylab("-log p-value")
  
vol

# logistic regression

gene_8 <- read.csv("C:\\Users\\Steven\\Desktop\\github\\LASSO\\8 gene_symbol.csv", head = T)
gene8_data <- subset(gene_data, select = as.character(gene_8$Symbol))
gene8_data <- gene8_data[1:38,]
gene8_data <- scale(gene8_data)
cache_state <- c(rep(0, 16), rep(1, 17), rep(0, 5))
gene8_data <- data.frame(gene8_data, cache_state)

model_logistic <- glm(cache_state ~ ., data = gene8_data, family = binomial(link = "logit"))
summary(model_logistic)

# LASSO

library(glmnet)

model_LASSO <- glmnet(x = as.matrix(gene8_data[,1:8]), y = gene8_data[,9], family = "binomial", alpha = 1)
plot(model_LASSO, xvar='lambda', main="LASSO", col = 1:8, lwd = 2)
legend("topright", lwd = 2, col = 1:8, legend = colnames(gene8_data[,1:8]), cex = .8)
abline(h = 0)

model_LASSO[["beta"]]

# cross validation
set.seed(10)
model_cv <- cv.glmnet(x = as.matrix(gene8_data[,1:8]), y = gene8_data[,9], family = "binomial", alpha = 1, nfolds = 10, type.measure = "mae")
plot(model_cv)
model_cv$lambda.min
log(model_cv$lambda.min)

plot(model_LASSO, xvar='lambda', main="LASSO", col = 1:8, lwd = 2)
legend("topright", lwd = 2, col = 1:8, legend = colnames(gene8_data[,1:8]), cex = .8)
abline(h = 0)
abline(v = log(model_cv$lambda.min), col = "red", lty = "dashed")

# post-LASSO

model_p <- glm(cache_state ~ S100A8 + S100A9 + S100A12 + CCL18, data = gene8_data, family = binomial(link = "logit"))
summary(model_p)

