% CREATES SIMULATED DATE FOR TRIAL SELECTION M3


load options



xi = options(:,1); % x = ref option
xj = options(:,2);
yi = options(:,3); % y = non-ref option
yj = options(:,4);

rx = (xi-xj)>=0; % this is for the Charness-Rabin model
zx = (xj-xi)>0;
ry = (yi-yj)>=0;
zy = (yj-yi)>0;

options(:,5) = rx; % x = ref option
options(:,6) = zx;
options(:,7) = ry;
options(:,8) = zy;


genmodels = [1 2 3]; % CHARNES-RABIN, COBB-DOUGLAS, RAWLSIAN
fitmodels = [1 2 3];

npars_vec = [3 2 2]; % number of parameters of each model

n_models = numel(genmodels);

model_agent_combination = cell(n_models,1); % simulated data is going to go here. 


R = 10; % number of repetitions per Unique trial
n_UT = size(options,1); % number of unique trials

tic
% create data structure
for m =1:n_models
    switch m
        case 1
            npars = npars_vec(m);
            pars_range = [1,10;... % alpha goes from 1 to 10 (noise)
                -1 1;... % weight to other when better off - ro in charness-rabin
                -1 1]; % weight to other when worse off - sigma in charness-rabinc
            alpha_vec = linspace(pars_range(1,1),pars_range(1,2),10); % FULL VERSION
            beta_vec = linspace(pars_range(2,1),pars_range(2,2),11);
            gamma_vec = linspace(pars_range(3,1),pars_range(3,2),11);
            
            model_agent_combination{m}.all_agents = combvec(alpha_vec,beta_vec,gamma_vec)';
            model_agent_combination{m}.n_agents = size(model_agent_combination{m}.all_agents,1);
            model_agent_combination{m}.agents = cell(model_agent_combination{m}.n_agents,1);
            
        case 2
            npars = npars_vec(m);
            pars_range = [1,10;...
                0 1];
            alpha_vec = linspace(pars_range(1,1),pars_range(1,2),35);
            beta_vec = linspace(pars_range(2,1),pars_range(2,2),35);
            
            model_agent_combination{m}.all_agents = combvec(alpha_vec,beta_vec)';
            model_agent_combination{m}.n_agents = size(model_agent_combination{m}.all_agents,1);
            model_agent_combination{m}.agents = cell(model_agent_combination{m}.n_agents,1);
            
        case 3
            npars = npars_vec(m);
            pars_range = [1,10;...
                0 1];
            alpha_vec = linspace(pars_range(1,1),pars_range(1,2),35);
            beta_vec = linspace(pars_range(2,1),pars_range(2,2),35);
            
            model_agent_combination{m}.all_agents = combvec(alpha_vec,beta_vec)';
            model_agent_combination{m}.n_agents = size(model_agent_combination{m}.all_agents,1);
            model_agent_combination{m}.agents = cell(model_agent_combination{m}.n_agents,1);
    end
end

% create DeltaU for each subject, each choice options
for m=1:n_models
    n_agents = model_agent_combination{m}.n_agents;
    for i = 1:n_agents
        pars = model_agent_combination{m}.all_agents(i,:);
        model_agent_combination{m}.agents{i}.pars = pars;
        switch m
            case 1 % CHARNESS-RABIN
                ro = pars(2); % How much I care about you when I'm better off
                sig = pars(3); % How much I care about you when I'm worse off
                
                DeltaU = (1 - ro.*rx - sig.*zx).*xi + (ro.*rx + sig.*zx).*xj - ((1 - ro.*ry - sig.*zy).*yi + (ro.*ry + sig.*zy).*yj);
                
            case 2 % COBB-DOUGLAS
                beta = pars(2); % How much I care about myself
                
                DeltaU = xi.^beta.*xj.^(1-beta) - (yi.^beta.*yj.^(1-beta));
                
            case 3 % RAWLSIAN
                beta = pars(2);
                
                DeltaU = min(beta*xi,xj) - min(beta*yi,yj);
        end
        
        
        model_agent_combination{m}.agents{i}.all_options = [options DeltaU];
    end
end

% create choice
for t= 1:n_UT
    for m=1:n_models
        n_agents = model_agent_combination{m}.n_agents;
        for i = 1:n_agents
            UT = model_agent_combination{m}.agents{i}.all_options(t,:);
            DeltaU = model_agent_combination{m}.agents{i}.all_options(t,9);
            alpha = model_agent_combination{m}.agents{i}.pars(1);
            py = 1./(1+exp(alpha*DeltaU));
            ct = rand(R,1) < py;
            ct = 2*ct - 1; % 1=Non-Ref-y; -1=Ref option-x
            model_agent_combination{m}.agents{i}.choice_data(:,:,t) = [repmat(UT,R,1) ct];
        end
    end
end


save('simulated_data_1.mat','model_agent_combination','-v7.3')


toc

% 25 min for the full version

