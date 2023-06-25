clear;
clc;
addpath(genpath('.'));
load('data/birds.mat');
    optmParameter.alpha   = 10^-3;                                                                                                                                                                                                                                                                              
    optmParameter.gamma   = 10^-3;   
    optmParameter.tuneParaOneTime = 1; 
    optmParameter.theta   = 100;
    optmParameter.maxIter            = 100;
    optmParameter.minimumLossMargin  = 10^-5;
    optmParameter.bQuiet             = 1;
    modelparameter.crossvalidation    = 1; 
    modelparameter.cv_num             = 5;
    modelparameter.L2Norm             = 1; 
    modelparameter.tuneParaOneTime    = 1;
    modelparameter.repetitions        = 1;
   starttime = datestr(now,0);

if exist('train_data','var')==1
    data=[train_data;test_data];
    target=[train_target,test_target];
    clear train_data test_data train_target test_target
end
data     = double(data);
target(target==-1)=0;
num_data = size(data,1);
if modelparameter.L2Norm == 1
    temp_data = normalization(data, 'l2', 1);
else
    temp_data = data;
end
clear data  
randorder = randperm(num_data);
Result_CNSF  = zeros(5,modelparameter.cv_num);

for i = 1:modelparameter.repetitions   
    for j = 1:modelparameter.cv_num
            fprintf('- Repetition - %d/%d,  Cross Validation - %d/%d\n', i, modelparameter.repetitions, j, modelparameter.cv_num);
            [cv_train_data,cv_train_target,cv_test_data,cv_test_target ] = generateCVSet( temp_data,target',randorder,j,modelparameter.cv_num );
            cv_train_target=cv_train_target';
            cv_test_target=cv_test_target';
            
            [model_CNSF]  = CNSF( cv_train_data, cv_train_target',optmParameter);
            Outputs=cv_test_data*model_CNSF;           
            Outputs=Outputs';
            Pre_Labels=sign( Outputs);
             Pre_Labels(Pre_Labels==-1)=0;
            Result_CNSF(:,j) = EvaluationAll(Pre_Labels,Outputs,cv_test_target);
    end
end
Avg_Result = zeros(5,2);
Avg_Result(:,1)=mean(Result_CNSF,2);
Avg_Result(:,2)=std(Result_CNSF,1,2);
fprintf('\nResults of CNSF\n');
PrintResults(Avg_Result);
rmpath(genpath('.'));
endtime = datestr(now,0);
% 