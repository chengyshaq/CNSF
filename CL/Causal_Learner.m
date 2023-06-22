function [Result1] = Causal_Learner(input_alg_name,data,data_type,alpha)

addpath(genpath(pwd));
maxK=3;
[samples,p]=size(data);
ns=max(data);
alg_name='PCstable';
total_alg=alg_name;    
if ~ismember(input_alg_name,total_alg)
    fprintf('\n%s is not a valid algorithm name \n',input_alg_name);
    return;
end

if strcmp(data_type,'dis')

    algorithm=str2func(strcat(input_alg_name,'_G2'));
    DAG=algorithm(data,alpha,ns,p,maxK);

elseif strcmp(data_type,'con')

    algorithm=str2func(strcat(input_alg_name,'_Z'));  
    DAG=algorithm(data,alpha,samples,p,maxK);   
end
    Result1=DAG;
end
    
    
    
    
    
    
    
    
    
    

