
function [proportion_chose_max] = proportion_chose_max(D)

xi = D(:,1);
xj = D(:,2);
yi = D(:,3);
yj = D(:,4);
ct = D(:,5);

X = xi+xj;
Y = yi+yj;

xmax = X > Y;
max_option = xmax*2-1;

chose_max = max_option==ct;

proportion_chose_max = sum(chose_max)/numel(chose_max);
        
        
