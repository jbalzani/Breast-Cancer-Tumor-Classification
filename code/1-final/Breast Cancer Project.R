options(digits = 3)
library(matrixStats)
library(tidyverse)
library(caret)
library(dslabs)
data("brca")

#q1
#proportion of malignant samples
mean(brca$y == "M")

#col number with highest mean
which.max(colMeans(brca$x))

#col with lowest std dev
sds <- rep(NA, ncol(brca$x))
for (i in 1:ncol(brca$x)) {
  sds[i] <- sd(brca$x[, i])
  }
which.min(sds)
which.min(colSds(brca$x))#better way

#q2 scale the matrix 
#use sweep 2 times to scale each col, subtract col mean, then divide by column's std dev
x_centered <- sweep(brca$x, 2, colMeans(brca$x), FUN = "-")
x_scaled <- sweep(x_centered, 2, colSds(brca$x), FUN = "/")
sd(x_scaled[, 1])
median(x_scaled[, 1])

#q3 calculate dist between all samples using the scaled matrix
d <- dist(x_scaled)
#avg dist between first sample (benign) and all other benign samples
d <- as.matrix(d)
#this works bc d is an nxn matrix, col #'s correspond to row #s
mean(d[1, brca$y == "B"])
mean(d[1, brca$y == "M"])

#official code
d_samples <- dist(x_scaled)
dist_BtoB <- as.matrix(d_samples)[1, brca$y == "B"]
mean(dist_BtoB[2:length(dist_BtoB)])
dist_BtoM <- as.matrix(d_samples)[1, brca$y == "M"]
mean(dist_BtoM)

#q4 make a heatmap of the relationship between features using scaled matrix
a <- t(x_scaled)
library(RColorBrewer)
heatmap(as.matrix(dist(a)),
        labRow = NA, labCol = NA)
#official code
d_features <- dist(t(x_scaled))
heatmap(as.matrix(d_features), labRow = NA, labCol = NA)

#q5 perform hierarchical clustering on the 30 features. cut tree into 5 groups
clust_algo <- hclust(d_features)
groups <- cutree(clust_algo, k = 5)
grp_tbl <- tibble(labels = clust_algo$labels, group = groups)

#official code
h <- hclust(d_features)
groups <- cutree(h, k = 5)
split(names(groups), groups)

#q6 perform pca of scaled matrix
pc <- prcomp(x_scaled)
var_explained <- cumsum(pc$sdev^2/sum(pc$sdev^2))

#q7 plot 1st 2 principal components with color being tumor type(benign/malig)
pc_df <- data.frame(pc$x[, 1:2], type = brca$y)
pc_df %>% ggplot(aes(x = PC1, y = PC2, color = type)) +
  geom_point()

#q8 make boxplot of 1st 10 PCs grouped by tumor type
for (i in 1:10) {
  boxplot(pc$x[, i] ~ brca$y, main = paste("PC", i))
}

#q9 create training and test sets
set.seed(1, sample.kind = "Rounding")  
test_index <- createDataPartition(brca$y, times = 1, p = 0.2, list = FALSE)
test_x <- x_scaled[test_index,]
test_y <- brca$y[test_index]
train_x <- x_scaled[-test_index,]
train_y <- brca$y[-test_index]

#check that test and training sets have similar proportion of benign tumors
mean(train_y == "B")
mean(test_y == "B")

#q10a perform k-means clustering on the training set with 2 centers & assign to k
#use predict_kmeans to make predictions on test set
set.seed(3, sample.kind = "Rounding")
predict_kmeans <- function(x, k) {
  centers <- k$centers    # extract cluster centers
  # calculate distance to cluster centers
  distances <- sapply(1:nrow(x), function(i){
    apply(centers, 1, function(y) dist(rbind(x[i,], y)))
  })
  max.col(-t(distances))  # select cluster with min distance to center
}
k <-  kmeans(train_x, centers = 2)
yhat_kmeans <- ifelse(predict_kmeans(test_x, k) == 2, "M", "B")
mean(yhat_kmeans == test_y)

#q10b what proportion of benign tumors are correctly id'ed?
#what proportion of malignant tumors are correctly id'ed?
sum(yhat_kmeans == "B" & test_y == "B")/sum(test_y == "B")
sum(yhat_kmeans == "M" & test_y == "M")/sum(test_y == "M")

#official code
sensitivity(factor(kmeans_preds), test_y, positive = "B")
sensitivity(factor(kmeans_preds), test_y, positive = "M")

#q11 fit logistic reg model on training set using all predictors. ignore warnings
#what is the accuracy on the test set?
train_data <- data.frame(train_x) %>%
  mutate(y = train_y)
reg_logit <- train(y ~ ., method = "glm", data = train_data)
yhat_logit <- predict(reg_logit, newdata = test_x)
mean(yhat_logit == test_y)

#q12 train a lda model and a qda model. what are the accuracies on the test set?
#lda model
lda_algo <- train(y ~ ., method = "lda", data = train_data)
yhat_lda <- predict(lda_algo, newdata = test_x)
mean(yhat_lda == test_y)

#qda model
qda_algo <- train(y ~ ., method = "qda", data = train_data)
yhat_qda <- predict(qda_algo, newdata = test_x)
mean(yhat_qda == test_y)

#q13 set seed to 5. fit loess model on training set. use default tuning grid.
#what is the accuracy on the test set?
set.seed(5, sample.kind = "Rounding")
modelLookup("gamLoess")
loess_algo <- train(y ~ ., method = "gamLoess", data = train_data)
yhat_loess <- predict(loess_algo, newdata = test_x)
mean(yhat_loess == test_y)

#q14 knn model. train knn model on training set using odd values of k from 3 to 21
#what is the final value of k used in the model?
#what is the accuracy of the knn model on the test set?
set.seed(7, sample.kind = "Rounding")
tuning <- data.frame(k = seq(3, 21, 2))
train_knn <- train(train_x, train_y,
                   method = "knn", 
                   tuneGrid = tuning)
train_knn$bestTune

knn_preds <- predict(train_knn, test_x)
mean(knn_preds == test_y)
##note: knn method in caret gives different answer than knn3() function.

#q15 train random forest model on training set using caret package.
#test mtry values of 3, 5, 7, 9
#use argument importance = TRUE
set.seed(9, sample.kind = "Rounding")
tuning <- data.frame(mtry = c(3, 5, 7, 9))
rf_algo <- train(train_y, train_x, method = "rf", importance = TRUE, 
                 tuneGrid = tuning)

#what value of mtry gives highest accuracy?
rf_algo$bestTune
#what is the accuracy on the test set?
yhat_rf <- predict(rf_algo, newdata = test_x)
mean(yhat_rf == test_y)
#what is the most important variable in the rf model?
varImp(rf_algo)

#q16 create an ensemble prediction to classify tumors as benign or malignant
#use k-means, logistic regression, lda, qda, loess, knn, random forest
#what is the accuracy of the ensemble prediction
ensemble_tbl <- tibble(yhat_kmeans = yhat_kmeans,
                        yhat_logit = yhat_logit,
                        yhat_lda = yhat_lda,
                        yhat_qda = yhat_qda,
                        yhat_loess = yhat_loess,
                        yhat_knn = knn_preds,
                        yhat_rf = yhat_rf) 

yhat_ensemble <- rep(NA, nrow(ensemble_tbl))
for (i in 1:nrow(ensemble_tbl)) {
  yhat_ensemble[i] <- ifelse(sum(ensemble_tbl[i,] == "M") >= 4, "M", "B")
}
mean(yhat_ensemble == test_y)

#official code
ensemble <- cbind(glm = glm_preds == "B", 
                  lda = lda_preds == "B", 
                  qda = qda_preds == "B", 
                  loess = loess_preds == "B", 
                  rf = rf_preds == "B", 
                  knn = knn_preds == "B", 
                  kmeans = kmeans_preds == "B")

ensemble_preds <- ifelse(rowMeans(ensemble) > 0.5, "B", "M")
mean(ensemble_preds == test_y)

#q16 make a table of the prediction of each model, and the ensemble prediction.
#what is most accurate?
#this does not work for some reason
ensemble_tbl <- ensemble_tbl %>%
  mutate(yhat_ensemble = as.factor(yhat_ensemble),
         yhat_kmeans = as.factor(yhat_kmeans))
accuracy <- rep(NA, ncol(ensemble_tbl))
for (i in 1:ncol(ensemble_tbl)) {
  accuracy[i] <- mean(ensemble_tbl[, i] == test_y)
}
#official code
models <- c("K means", "Logistic regression", "LDA", "QDA", "Loess", "K nearest neighbors", "Random forest", "Ensemble")
accuracy <- c(mean(kmeans_preds == test_y),
              mean(glm_preds == test_y),
              mean(lda_preds == test_y),
              mean(qda_preds == test_y),
              mean(loess_preds == test_y),
              mean(knn_preds == test_y),
              mean(rf_preds == test_y),
              mean(ensemble_preds == test_y))
data.frame(Model = models, Accuracy = accuracy)
colnames(ensemble_tbl[, which.max(accuracy)])