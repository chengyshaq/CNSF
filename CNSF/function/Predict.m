function [TY,results] = Predict(modelTrain, Xt, Yt)
TY = Xt*modelTrain;
results =  evalt(TY', Yt, (max(TY(:))-min(TY(:)))/2, 1);
end    