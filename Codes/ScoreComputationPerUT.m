function ScoreComputationPerUT(j)

tic
load simulated_data_1.mat
toc
n_UT = size(model_agent_combination{1}.agents{1}.all_options,1);

n_models = 3;
R = size(model_agent_combination{1}.agents{1}.choice_data,1); % Repetitions per Unique Trial


% We sent batches of 50 trials to the cluster (High Performance Computer)
% to calculate the mutual information between models (>1000 agents each)
% and their responses. Each batch corresponds to a "j" value, which go from
% 1 to 59

it = (j-1)*50 +1; 
ft = j*50;
if ft == 3000;
    ft = 2934;
end

% A version of this below went to the cluster, for each of the 59 batches. 

SCORE = nan(ft-it+1,1);

tic
% loop through unique trials
for tidx = 1:ft-it+1
    t = it-1+tidx;
    % initialize:
    SCORE_model_agent = cell(n_models,1); % logposterior of model given response for EACH agent 
    SCORE_model = nan(1,n_models); % logposterior of each model given responses to this trial accross ALL agents of each model
    
    for m=1:n_models
        n_agents = model_agent_combination{m}.n_agents;
        SCORE_model_agent{m} = nan(1,n_agents); 
        for i = 1:n_agents
            ct = model_agent_combination{m}.agents{i}.choice_data(:,10,t); % column 10 contain choice data, ct. 
            % initialize:
            likelihood_model_agent = cell(n_models,1);
            likelihood_models = nan(R,n_models); % marginalizing over agents within model. R repetitions of the same trial
            for hm=1:n_models % hypothesized models
                n_h_agents = model_agent_combination{hm}.n_agents;
                likelihood_model_agent{hm} = nan(R,n_h_agents); % likelihood of hypothesized model by agent given response
                for hi = 1:n_h_agents % hypothesized agents
                    DeltaU = model_agent_combination{hm}.agents{hi}.all_options(t,9); % colum 9 contains Delta Utility between options on this trial, for this agent
                    alpha = model_agent_combination{hm}.agents{hi}.pars(1); % alpha here is the noise parameter
                    likelihood_model_agent{hm}(:,hi) = 1./(1+exp(ct*alpha*DeltaU));
                end
                likelihood_models(:,hm) = sum(likelihood_model_agent{hm},2)*(1/n_h_agents); % marginalizing over agents to get likelihood of model ; assuming uniform distribution of agents within model
            end
            % log posterior probability of true model
            log_posterior_current_model = log(likelihood_models(:,m)./sum(likelihood_models,2)); % p(M|r) = p(r|M)*p(M)/p(r); then sum 
            % average log posterior prob of current model over all repetitions 
            SCORE_model_agent{m}(i) = mean(log_posterior_current_model);
            
        end
        % average scores over all agents of same model
        SCORE_model(m) = mean(SCORE_model_agent{m});
    end
    
    SCORE(tidx) = mean(SCORE_model, 2);
    %save test.mat
end
toc

filename = strcat('Score_',num2str(j)); 
save(filename, 'SCORE')
