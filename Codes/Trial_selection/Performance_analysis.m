% PERFORMANCE

% evaluate performance with 100 selected trials.
% then evaluate 80, 200, 300, and 500 selected trials. 

%%% A version of this was sent to the Cluster (High Performance Computer)
%%% to carry out 10 repetitions of this performance analysis for each number
%%% of selected trials. 


load SCORES.mat 

% take the data
load simulated_data_3.mat
n_UT = size(model_agent_combination{1}.agents{1}.all_options,1);

n_models = 3;
n_selected_trials = 80; % change this number to [50 80 100 300 500];
selected_trials = IX(1:n_selected_trials);

% randomely take one observation of each trial
tic
new_data = cell(n_models,1);
for m=1:n_models
    n_agents = model_agent_combination{m}.n_agents;
    new_data{m}.n_agents = n_agents;
    new_data{m}.agents = cell(n_agents,1);
    new_data{m}.all_agents = model_agent_combination{m}.all_agents;
    for i = 1:n_agents
        new_data{m}.agents{i}.choice_data = nan(n_selected_trials,10); 
        for tidx=1:n_selected_trials
            t = selected_trials(tidx); % t is the true trial index in the original pool of ~3000 trials.
            new_data{m}.agents{i}.choice_data(tidx,:) = model_agent_combination{m}.agents{i}.choice_data(randsample(10,1),:,t); % randonmly take 1 of the R repetitions of the same trial
        end
    end
end
toc 

N_AGENTS = nan(1,n_models); % n_agents of the 3 models
for i=1:n_models
    N_AGENTS(i) = model_agent_combination{i}.n_agents;
end

WINNING_MODELS = nan(n_models,max(N_AGENTS));


models_names = {'FEHR-SCHMIDT','UTILITARIAN','RAWLSIAN','COBB-DOUGLAS','CHARNES-RABIN','CES','LEXICOGRAPHIC','COX_CES'};
fitmodels = [5 4 3]; % these later changed to 1,2 and 3 when we decided to focus on Rawlsian, Cobb-Douglas and Charness-Rabin models. 1,2,3 notation appears in "ealier" codes: codes used earlier within this trial selection process
nfmodels = length(fitmodels);
maxnumpars = 4;
nparshat_vec = [3 2 2 2 3 3 1 4]; 

initialVars = structure_initial_variables (fitmodels, nfmodels, maxnumpars, nparshat_vec); 

% create fit results structure:
tic
new_data_Fit = cell(n_models,1);
for m=1:n_models
    n_agents = model_agent_combination{m}.n_agents;
    new_data_Fit{m}.n_agents = n_agents;
    new_data_Fit{m}.agents = cell(n_agents,1);
    for i = 1:n_agents
        new_data_Fit{m}.agents{i}.AIC = nan(1,n_models);
        new_data_Fit{m}.agents{i}.parshat = nan(1,maxnumpars);
        S_DATA = new_data{m}.agents{i}.choice_data;
        D = [S_DATA(:,1:4) S_DATA(:,10) S_DATA(:,5:8)]; % results in: [xi xj yi yj ct rx zx ry zy]
        [AIC winning_model winning_params] = Fit_sim_subject (D, initialVars); 
        new_data_Fit{m}.agents{i}.AIC = AIC;
        new_data_Fit{m}.agents{i}.parshat = winning_params;
        WINNING_MODELS(m,i) = winning_model;
    end
end
toc
% 1,15 hs for 100 trials 1 repetition - Data is saved as: new_sim_data_1repetition_100trials.mat
% 1,5 hs for 300 trials 1 repetition
% 2,12 hs for 500 trials 1 repetition


p1 = sum(WINNING_MODELS(1,1:N_AGENTS(1))==5)/N_AGENTS(1); % model 5 = CHARNES-RABIN
p2 = sum(WINNING_MODELS(2,1:N_AGENTS(2))==4)/N_AGENTS(2); % model 4 = COBB-DOUGLAS
p3 = sum(WINNING_MODELS(3,1:N_AGENTS(3))==3)/N_AGENTS(3); % model 3 = RAWLSIAN

[p1;p2;p3]


% 100 trials 1 repetition:
%0.9083 CH-R
%0.9820 C-D
%0.9192 R
% 300 trials 1 repetition:
%0.9195
%0.9580
%0.9371
% 500 trials 1 repetition:
% 0.9281
% 0.9910
% 0.9616

% CHANGE TO NEW NOMENCLATURE OF MODELS:
WINNING_MODELS(WINNING_MODELS(:) == 5) = 1;
WINNING_MODELS(WINNING_MODELS(:) == 4) = 2;
WINNING_MODELS(WINNING_MODELS(:) == 3) = 3;
% MODEL 1 = CHARNESS-RABIN
% MODEL 2 = COBB-DOUGLAS
% MODEL 3 = RAWLSIAN


% extract parameters --> COMPARE PARS TO PARSHAT

npars_vec = [3 2 2 ]; % BEWARE: this changes. Before it included 8 models
nparshat_vec = [3 2 2 ];
maxnumpars = 3; % BEWARE: this changes. Before maximum was 4. 

PARS = nan(sum(N_AGENTS),maxnumpars+1); % add a column for the true model;
PARSHAT = nan(sum(N_AGENTS),maxnumpars+1); % add a column fot the winning model;

j=0;
for m=1:n_models
    n_agents = N_AGENTS(m);
    for i=1:n_agents
        j=j+1;
        PARS(j,1:npars_vec(m)) = new_data{m}.all_agents(i,:); 
        PARS(j,maxnumpars+1) = m; % adds column of true model
        winning_model = WINNING_MODELS(m,i);
        PARSHAT(j,1:nparshat_vec(winning_model)) = new_data_Fit{m}.agents{i}.parshat(1:nparshat_vec(winning_model)); % this was necessary cause of error in new_data_fit maxnumpar. used old 4 instead of new 3. 
        PARSHAT(j,maxnumpars+1) = winning_model;
    end
end


genmodels = [1 2 3];
n_genmodels = numel(genmodels);
for g = 1:n_genmodels
    npars = npars_vec(g);
    gs = find([PARS(:,maxnumpars+1) == g].*[PARSHAT(:,maxnumpars+1) == g]);
    
    spars = PARS(gs, 1:npars);
    sparshat = PARSHAT(gs,1:npars);
    
    figure(g)
    for p = 1:npars
        subplot(1,npars,p)
        scatter(spars(:,p),sparshat(:,p)); xlabel('true par'); ylabel('estimated par'); hold on
        plot([min(spars(:,p)) max(spars(:,p))], [min(spars(:,p)) max(spars(:,p))],'k--')
        suptitle(sprintf('subjects in categorty %d',g))
    end
end
hold off


filename = sprintf('new_sim_data_1repetition_%dtrials',n_selected_trials);
save(filename,'new_data','new_data_Fit','p1','p2','p3','WINNING_MODELS','PARSHAT')


%% Compile results from the Cluster

N_selected_trials = [50 80 100 300 500];

for k=1:length(N_selected_trials)
    n_selected_trials = N_selected_trials(k);
    PERFORMANCE = nan(3,10); % 10 repetitions of the performance analysis for each N_selected_trials. Reported results are an average of this. 

    folder = fullfile('Performance',[sprintf('%d',n_selected_trials),'_trials']);

    for j=1:10
        filename = sprintf('performance_%d.mat',j);
        load(fullfile(folder,filename))
        PERFORMANCE(:,j)=performance;
    end

    mean_performance = mean(PERFORMANCE,2);

    save(sprintf('PERFORMANCE_%d_trials',n_selected_trials),'PERFORMANCE','mean_performance')

end

%% GET PERFORMANCE WITH RANDOM WLARASIAN TRIALS (downward facing trials)

load simulated_data_3.mat
load options.mat

DS = options(:,3) - options(:,1);
DO = options(:,4) - options(:,2);
TI = abs(DS)>0.01 & abs(DO)>0.01;
TII = DO./DS<0;
TIII = DO./DS~=Inf;
TIV = TI & TII & TIII;
w_options = options(TIV,:); % walrasian options = downward facing trials (positive price, positive sacrifice for the other to get something)

figure
scatter(options(:,4),options(:,3),10,'oc','filled');hold on
scatter(options(:,2),options(:,1),10,'ob','filled');
for t=1:sum(TIV)
    plot([w_options(t,2),w_options(t,4)],[w_options(t,1),w_options(t,3)],'-m')
    scatter(w_options(t,4),w_options(t,3),15,'om')
end

hola = find(TIV); % get 10 random samples for the 10 repetitions of the perfomrance analysis. 
for i=1:10
    r_idx(:,i) = randsample(hola,25); % we chose 25 trials to approximate the number of trials typically used in other experiments (although they may repeat them)
end


save('r_idx','r_idx')

% Go to the Performance_analysis code and define selected_trials=r_idx(:,j)
% with j=1:10 for the 10 repetitions.


% COMPILE RESULTS
n_models = 3;
n_jobs = 10;
folder = '25_rand';
PERFORMANCE = nan(n_models, n_jobs);
for j=1:10
    filename = strcat('performance_',num2str(j));
    load(fullfile(folder,filename));
    PERFORMANCE(:,j) = performance;
end

mean_performance = mean(PERFORMANCE,2);
save(strcat('PERFORMANCE_',folder), 'PERFORMANCE','mean_performance')
    
%% NOW REPEAT BUT USING AN EGALITARIAN REFERENCE POINT. 45,45

TV = options(:,1)==45 & options(:,2)==45;
TVI = TI & TII & TIII & TV;

r_idx_4545 = find(TVI);
r_idx_4545(r_idx_4545>2304) = [];

save('r_idx_4545','r_idx_4545')

% COMPILE RESULTS
n_models = 3;
n_jobs = 10;
folder = '25_rand_4545';
PERFORMANCE = nan(n_models, n_jobs);
for j=1:10
    filename = strcat('performance_',num2str(j));
    load(fullfile(folder,filename));
    PERFORMANCE(:,j) = performance;
end

mean_performance = mean(PERFORMANCE,2);
save(strcat('PERFORMANCE_',folder), 'PERFORMANCE','mean_performance')
 

