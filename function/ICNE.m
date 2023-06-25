function [LM] = ICNE( train_target )

[n,nc]=size(train_target);
for i = 1:nc
    sum1 = sum(train_target(:,i)==1);
    sum2 = sum(train_target(:,i)==0);
    P =   sum1/n + 1;
    q =   -sum2/n  ;
    
    for j = 1:n
        if train_target(j,i) == 1
           train_target(j,i)=P;
        else
          train_target(j,i)=q;         
        end
    end  
end
LM = train_target;

end