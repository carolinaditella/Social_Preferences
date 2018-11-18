function error = loglike_collusion (pars, D)


xi = D(:,1);
xj = D(:,2);
yi = D(:,3);
yj = D(:,4);
ct = D(:,5);

rx = D(:,6);
zx = D(:,7);
ry = D(:,8);
zy = D(:,9);

a = pars(1);


% LINEAR 
w = 1;

DeltaU = xi + w*xj - (yi + w*yj);
LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));

    
error = - LOGS;

end

