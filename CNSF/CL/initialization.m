

function [optmParameter,modelparameter] =  initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Optimization Parameters
    optmParameter.alpha   = 10^-15;                                                                                                                                                                                                                                                                               ; % 2.^[-10:10] % sparsity
    optmParameter.beta    = 10^-8;  
    optmParameter.gamma      = 10^-8;      

    optmParameter.searchPara      = 0; % indicate whether tuning the parameters, {0:not,1:yes}
    optmParameter.tuneParaOneTime = 1; % indicate that tuning the parameter one time or tuning it in each fold. {0: each fold,1: only one time}
    optmParameter.theta   = 10;

    optmParameter.maxIter            = 100;
    optmParameter.minimumLossMargin  = 10^-5;
    optmParameter.bQuiet             = 1;

    %% Model Parameters
    modelparameter.crossvalidation    = 1; % {0,1}
    modelparameter.cv_num             = 5;
    modelparameter.L2Norm             = 1; % {0,1}

    modelparameter.tuneParaOneTime    = 1;
    modelparameter.repetitions        = 1;

end