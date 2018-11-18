function error = loglike (pars, D, fitmodel)


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


switch fitmodel
    case 1 % FEHR-SCHMIDT
        d = pars(2); % envy
        b = pars(3); % guilt
        
        DeltaU = xi - d*max(xj-xi,0) - b*max(xi-xj,0) - (yi - d*max(yj-yi,0) - b*max(yi-yj,0) );
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
    case 2 % LINEAR (UTILITARIAN)
        w = pars(2);
        
        DeltaU = xi + w*xj - (yi + w*yj);
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
        
    case 3 % RAWLSIAN
        w = pars(2);
        
        DeltaU = min(w*xi,xj) - min(w*yi,yj);
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
    case 4 % COBB-DOUGLAS
        w = pars(2);
        
        DeltaU = xi.^w.*xj.^(1-w) - (yi.^w.*yj.^(1-w));
        LOGS = sum(-log(1 + exp(ct*a.*DeltaU)));
        
    case 5 % CHARNESS-RABIN
        ro = pars(2);
        sig = pars(3);
        
        DeltaU = (1 - ro.*rx - sig.*zx).*xi + (ro.*rx + sig.*zx).*xj - ((1 - ro.*ry - sig.*zy).*yi + (ro.*ry + sig.*zy).*yj);
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
    case 6 % CES
        w = pars(2);
        ro = pars(3);
        
        DeltaU = ((1-w)*xi.^ro + w*xj.^ro).^(1/ro) - ((1-w)*yi.^ro + w*yj.^ro).^(1/ro);
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
    case 7 % LEXSELF
        
        DeltaU = xi + xj./(xj + 1) -  (yi + yj./(yj + 1));
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
    case 8 % COX_CES
        tu = pars(2);
        td = pars(3);
        ro = pars(4);
        
        DeltaU = (rx.*((1-tu)*xi.^ro + tu*xj.^ro).^(1/ro) + zx.*((1-td)*xi.^ro + td*xj.^ro).^(1/ro)) - (ry.*((1-tu)*yi.^ro + tu*yj.^ro).^(1/ro) + zy.*((1-td)*yi.^ro + td*yj.^ro).^(1/ro));
        LOGS = sum(-log(1+ exp(ct*a.*DeltaU)));
        
end

error = - LOGS;
