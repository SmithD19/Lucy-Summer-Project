
sample <- sample.int(n = nrow(foo), size = floor(.80 * nrow(foo)), replace = F)
train <- foo[sample,]
test <- foo[-sample,]

library(xgboost)

sparse_mat <- sparse.model.matrix(west_nile_virus ~., -1, data = train)

output_vector = train$west_nile_virus

output_vector = train[,"west_nile_virus"] == 1

bst <- xgboost(data = data.matrix(train), label = output_vector, max.depth = 4,
               eta = 1, nthread = 4, nrounds = 100, objective = "binary:logistic", verbose = 1)

pred <- predict(bst, data.matrix(test))

prediction <- as.numeric(pred > 0.5)
print(head(prediction))

err <- mean(as.numeric(pred > 0.5) != test$west_nile_virus)
print(paste("test-error=", err))

#-----

dtrain <- xgb.DMatrix(data = data.matrix(train), label = train$west_nile_virus)
dtest <- xgb.DMatrix(data = data.matrix(test), label = test$west_nile_virus)

watchlist <- list(train=dtrain, test=dtest)

bst <- xgb.train(data = dtrain, max.depth = 2, eta = 1, nthread = 2, nrounds = 4, 
                 watchlist=watchlist, objective = "binary:logistic")


pred <- predict(bst, data.matrix(foo))

prediction <- as.numeric(pred > 0.5)
print(head(prediction))

err <- mean(as.numeric(pred > 0.5) != foo$west_nile_virus)
print(paste("test-error=", err))

xgb.importance(data = data.matrix(train), label = train$west_nile_virus, 
               model = bst)

xgb.importance(colnames(train), model = bst)



