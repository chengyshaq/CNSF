function [CL] = Causal(X)
    data_type='dis';
    alg_name = 'PCstable';
    data = X';
    data(data==1)=2;
    data(data==0)=1;
    CL = Causal_Learner(alg_name,data,data_type,0.05);
end