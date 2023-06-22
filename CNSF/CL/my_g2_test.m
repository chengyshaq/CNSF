function [pvalue, dep] = my_g2_test(X, Y, S, Data, ns, alpha)
[nSamples,nVars]=size(Data);
if isreliablecit(ns, nSamples, X, Y, S)
    [p,stat] = citpvalue(Data,ns, X, Y, S);  
    if p > alpha        
        ca = 0;
        ca2 = 0;       
    else        
        ca = 1/p;
        ca2 = abs(stat);      
    end
else  
    p = NaN;
    stat = NaN;
    if isempty(S)
        ca = 0;
        ca2 = 0;
    else
        ca = Inf;
        ca2 = Inf;
    end
end


pvalue=p;
dep=ca2;

end

function [p,stat] = citpvalue(Data, ns, i, j, k)
ind = [i j k];
testSample = Data(:, ind);
testNLevels = ns(:, ind);
nObs = size(testSample, 1);
nTestVars = size(testSample, 2);
Obs = accumarray(testSample, ones(1, nObs), testNLevels);
ObsSum2 = sum(Obs, 2);
Obs_xs_size = ones(1, nTestVars);
Obs_xs_size(2) = testNLevels(2);
Obs_xs = repmat(ObsSum2, Obs_xs_size); 

ObsSum1 = sum(Obs, 1); 
Obs_ys_size = ones(1, nTestVars);
Obs_ys_size(1) = testNLevels(1);
Obs_ys = repmat(ObsSum1, Obs_ys_size); 

Obs_s = sum(ObsSum1, 2); 
Obs_s_size = ones(1, nTestVars);
Obs_s_size([1 2]) = testNLevels([1 2]);
Obs_s = repmat(Obs_s, Obs_s_size);

Exp = Obs_xs.*Obs_ys./Obs_s;
j3NLevelsProd = prod(testNLevels(3:end));

ObsSum1 = reshape(ObsSum1, testNLevels(2), j3NLevelsProd);
ObsSum2 = reshape(ObsSum2, testNLevels(1), j3NLevelsProd);
df = 0;

for iComb = 1:j3NLevelsProd
    
    df = df + max(testNLevels(1) - 1 - sum(~ObsSum2(:,iComb)), 0) * max(testNLevels(2) - 1 - sum(~ObsSum1(:,iComb)), 0);
    
end

if df == 0
    p = 1;
    stat = 0;   
    return;    
end

Obs_vector = Obs(:);
Exp_vector = Exp(:);
stat = chi2stat(Obs_vector, Exp_vector);
p = gammainc(stat/2, df/2, 'upper'); 
end


function stat = chi2stat(obs, exp)
terms = obs.*log(obs./exp);
terms(isnan(terms)) = 0;
stat = 2*sum(terms);
end

function tf = isreliablecit(ns, nSamples,i, j, k)
hps=5; 
testVarNValues = ns([i j k]);
tf = nSamples / prod(testVarNValues) >= hps;
end
