% ANALYSIS CODE


% - FIT ALL SUBJECTS - Including tests for collusion and rationaity

% This has been adapted to run all subjects (Study 1 and 2) at the same
% time. 

define_initialVars

% define Analysis Variables:

%s_n_list = 301:330; % for study 1
% s_n_list = 331:360; % for study 2
s_n_list = 301:360; % for study 1 and 2. 

partners = {'friend','stranger'};
n_partners = length(partners);
n_subjects = length(s_n_list);

analysisVars = structure_analysis_variables(partners,n_partners,s_n_list,n_subjects); 


AS = cell(n_subjects,1); % AS = Analyzed Subjects. 

% CapsLocks typically refer to matrices with the info of every subject. 
WINNING_MODELS = nan(n_subjects,n_partners);
WINNING_PARAMETERS = nan(n_partners,maxnumpars,n_subjects);
TEMPORAL_AEI = nan(n_subjects,12); % 12 for maximum 12 blocks (in temporal order); regardless of partner
COLLUSION_info = nan(n_subjects,1); % AIC for params 0.5,0.5
COLLUSION_info_II = nan(n_subjects,1); % proportion of trials where chosen option is greater total sum


for i=1:length(s_n_list)
    
    s_n = s_n_list(i);
    
    data_folder = '/Users/carolinaditella/Dropbox/Social Preferences/Option Presentation/SOC_PREF_DATA'; % potentially change Option Presentation folder to DATA folder.
    datadir =  fullfile(data_folder,num2str(s_n),['SOC_PREF_DATA_', num2str(s_n),'.mat']);
    load(datadir) % Data (Raw)
    
   
    for j=1:n_partners
        partner = partners{j}
        
        [AICtemp LLtemp parshatMtemp winning_modeltemp winning_paramstemp Dtemp] = Fit_this_subject(Data, initialVars, analysisVars, partner);
        
        AIC(1,:,j) = AICtemp;
        LL(1,:,j) = LLtemp;
        parshatM(:,:,j) = parshatMtemp;
        winning_model(j) = winning_modeltemp;
        winning_params (j,:) = winning_paramstemp;
        D{j} = Dtemp;
        
        if j==1
            COLLUSION_info(i) = test_collusion_ga(Dtemp);
            COLLUSION_info_II(i) = proportion_chose_max(Dtemp);
        end
        
    end
    
    
    % RATIONALITY
    % this is outside of the partners' loop because it takes all blocks
    % -friend's and stranger's- in cronological order and then separates
    % them by partner. it remains within the subject's loop (i). 
    
    [mean_AEI, block_wise_AEI, temporal_AEI] = get_GARP(s_n);
    
    Data.rationality.block_wise_AEI = block_wise_AEI;
    Data.rationality.mean_AEI = mean_AEI;
    Data.rationality.temporal_AEI = temporal_AEI;
    
    % start putting every subject's AEI together:
    TEMPORAL_AEI(i,1:length(temporal_AEI)) = temporal_AEI;
    
    %%%%
    
    % STORE:
    Data.AIC = AIC; % this has both partners' info
    Data.LL = LL;
    Data.parshatM = parshatM;
    Data.winning_model = winning_model;
    Data.winning_params = winning_params;
    Data.D = D;
    
    AS{i} = Data;
    
    % Start putting every subject's winning model and winning parameters
    % together
    WINNING_MODELS(i,:) = winning_model;
    WINNING_PARAMETERS(:,:,i) = winning_params;
    
    
end


%save ('AS_CESS', 'AS','WINNING_MODELS','WINNING_PARAMETERS','COLLUSION_info','COLLUSION_info_II','TEMPORAL_AEI')
%save('AS_CESS.mat','initialVars','analysisVars,'-append')
%%

% PUT ALL THE AIC SCORES TOGETHER
AIC_MAT = nan(n_partners,nfmodels,n_subjects);
for i = 1:length(s_n_list)
    AIC_MAT(1,:,i) = AS{i}.AIC(1,:,1);
    AIC_MAT(2,:,i) = AS{i}.AIC(1,:,2);
end
AIC_F = squeeze(AIC_MAT(1,:,:))';
AIC_S = squeeze(AIC_MAT(2,:,:))';


AIC_COMBINED = AIC_F + AIC_S;


study_1=1:30; study_2=31:60; all_studies = 1:60;

% Compare models - number of subjects
sum(AIC_COMBINED(study_1,1)<AIC_COMBINED(study_1,2)) % compare Ch-R versus Cobb-Douglas. 24 gana CH-R
sum(AIC_COMBINED(study_1,1)<AIC_COMBINED(study_1,3)) % compare Ch-R versus Rawlsian. 26 gana CH-R
sum(AIC_COMBINED(study_1,2)<AIC_COMBINED(study_1,3)) % compare Cobb-Douglas versus Rawlsian. 23 gana C-D

% number of subject for which CH-R wins over both C-D and R
sum(AIC_COMBINED(study_1,1)<AIC_COMBINED(study_1,2) & AIC_COMBINED(study_1,1)<AIC_COMBINED(study_1,3)) %22

% for those subjects for whch CH-R looses against C-D and R, by how much. 
AIC_COMBINED(AIC_COMBINED(study_1,2)<AIC_COMBINED(study_1,1),1) - AIC_COMBINED(AIC_COMBINED(study_1,2)<AIC_COMBINED(study_1,1),2)
AIC_COMBINED(AIC_COMBINED(study_1,3)<AIC_COMBINED(study_1,1),1) - AIC_COMBINED(AIC_COMBINED(study_1,3)<AIC_COMBINED(study_1,1),3)

% diff AIC: 1- CH-R vs C-D, 2- CH-R vs R
dif = [AIC_COMBINED(study_1,2)- AIC_COMBINED(study_1,1), AIC_COMBINED(study_1,3)- AIC_COMBINED(study_1,1),AIC_COMBINED(study_1,3)- AIC_COMBINED(study_1,2)];
sum_diff = sum(dif);

k=30; % with boostrap confidence intervals. 
for j=1:1000
    sample_j = datasample(dif,k);
    sums(j,:) = sum(sample_j);
end

mean(sums);
CI = prctile(sums, [2.5 97.5])   

% strangers only
sum(AIC_S(study_1,1)<AIC_S(study_1,2))
sum(AIC_S(study_1,1)<AIC_S(study_1,3))
sum(AIC_S(study_1,2)<AIC_S(study_1,3))

% Compare CH-R with friends with collusion
sum(AIC_F(all_studies,1)<COLLUSION_info(all_studies))

% GROUP ANALYSIS FOR STUDY 1
evidence = 2*repmat([3,2,2],numel(study_1),1) - AIC_COMBINED(study_1,:)/2; % converts AIC to Loglikelihood again. 
[alpha,exp_r,xp,g,ll] = spm_BMS(evidence); 


%% FIT AND COMPARE THE SIMPLE LINEAR MODEL. 

AIC_simple_linear = nan(n_subjects,n_partners);
parshat_simple = nan(n_subjects,n_partners);

initialVars_simple = initialVars;
initialVars_simple.fitmodels = 2; % model number 2 is the simple linear
initialVars_simple.nfmodels = 1;

for s = 1:n_subjects
    
    s_n = s_n_list(s);
    
    data_folder = '/Users/carolinaditella/Dropbox/Social Preferences/Option Presentation/SOC_PREF_DATA';
    datadir =  fullfile(data_folder,num2str(s_n),['SOC_PREF_DATA_', num2str(s_n),'.mat']);
    load(datadir)
    
    for p=1:n_partners
        partner = partners{p}
        
        [AICtemp, ~, parshatMtemp, winning_modeltemp, ~, ~] = Fit_this_subject(Data, initialVars_simple,analysisVars, partner);
        
        AIC_simple_linear(s,p) = AICtemp;
        parshat_simple(s,p) = parshatMtemp(2);
        
    end
end
 
AIC_diff_simple_linear = sum(AIC_simple_linear,2) - AIC_COMBINED(:,1); % minus Charnes Rabin. Combined because it combines friend and stranger. In simple_linear it's done summing the columns. 

% - Only Study 1:
sum(AIC_diff_simple_linear(study_1)>0) % number of subjects for which CH-R is better than simple linear

sum(AIC_diff_simple_linear(study_1)) % by how much in total

% bootstrap confidence intervals for total sum
k=30
for j=1:1000
    sample_j = datasample(AIC_diff_simple_linear(study_1),k);
    SUMS(1,j) = sum(sample_j);
end

mean(SUMS) % 317.93 % remarkably similar to actual sum

CI = prctile(SUMS,[2.5 97.5]) % 47.90 747.50

% for strangers only
AIC_diff_simple_linear_strangers_only = AIC_simple_linear(:,2) - AIC_S(:,1); % minus Charnes Rabin. 

sum(AIC_diff_simple_linear_strangers_only(study_2)>0) % use study 2 for this comparison because it has better parameter estimates. 




%% - PSYCHOMETRIC DATA 

% Beware: the psychometric has 5 price levels in study 1 and 7 levels in study 2.

n_subjects = numel(all_studies);

PSYCHOMETRIC_DATA = cell(n_subjects,1);
for s=1:n_subjects
    PSYCHOMETRIC_DATA{s}.partner = Get_Psychometric(s,initialVars,analysisVars);
end

% By Status
for s=1:n_subjects
    for p=1:n_partners 
        PSYCHOMETRIC_DATA{s}.partner{p}.status = Get_Psychometric_status(s,p,analysisVars);
    end
end

%save('AS_CESS','PSYCHOMETRIC_DATA','-append')


% DONATION EFFICIENCY PRYCHOMETRIC (donation efficiency = inverse price of giving)

% plot for subjects in study 1: 
n_subjects = numel(study_1);
n_levels = 5;
n_models = nfmodels;


% values for all considered subjects:
INVERSE_PRICE_MEAN = nan(1,n_levels);
GEN_CHOICE = nan(n_subjects,n_levels,n_partners);
GEN_CHOICE_MEAN = nan(1,n_levels,n_partners);
GEN_CHOICE_SEM = nan(1,n_levels,n_partners);
MODEL_PREDICTION = nan(n_subjects,n_levels,n_models,n_partners);
MODEL_PREDICTION_MEAN = nan(1,n_levels,n_models,n_partners);
MODEL_PREDICTION_SEM = nan(1,n_levels,n_models,n_partners);
    
for s=1:n_subjects
    for p=1:n_partners        
        INVERSE_PRICE_MEAN = PSYCHOMETRIC_DATA{s}.partner{p}.Inverse_price_of_giving.collapsed.inverse_price_mean;
        GEN_CHOICE(s,:,p) = PSYCHOMETRIC_DATA{s}.partner{p}.Inverse_price_of_giving.collapsed.gen_choice_mean;
        for m=1:n_models
            MODEL_PREDICTION(s,:,m,p) = PSYCHOMETRIC_DATA{s}.partner{p}.Inverse_price_of_giving.collapsed.model_prediction_mean(m,:);
            MODEL_PREDICTION_MEAN(1,:,m,p) = mean(MODEL_PREDICTION(:,:,m,p),1);
            MODEL_PREDICTION_SEM(1,:,m,p) = std(MODEL_PREDICTION(:,:,m,p))./sqrt(size(MODEL_PREDICTION,1));
        end
    GEN_CHOICE_MEAN(1,:,p) = mean(GEN_CHOICE(:,:,p));
    GEN_CHOICE_SEM(1,:,p) = std(GEN_CHOICE(:,:,p))./sqrt(size(GEN_CHOICE,1));
    end
end

% PLOT IT:
C = [0 0 1; 0 1 0; 1 0 0];
partners={'FRIEND','STRANGER'};
tick_labels = {num2str(INVERSE_PRICE_MEAN(1),2),num2str(INVERSE_PRICE_MEAN(2),2),num2str(INVERSE_PRICE_MEAN(3),2),num2str(INVERSE_PRICE_MEAN(4),2),'Free'};
for p=1:2
    figure
    for m=1:n_models
        hf(m) = fill ([INVERSE_PRICE_MEAN'; flipud(INVERSE_PRICE_MEAN')],[MODEL_PREDICTION_MEAN(:,:,m,p)'- MODEL_PREDICTION_SEM(:,:,m,p)'; flipud(MODEL_PREDICTION_MEAN(:,:,m,p)'+ MODEL_PREDICTION_SEM(:,:,m,p)')],C(m,:),'LineStyle','none','FaceAlpha',0.3); hold on; %,'FaceAlpha',0.3 in case I want to do it transparent
    end
    hp = errorbar(INVERSE_PRICE_MEAN,GEN_CHOICE_MEAN(:,:,p),GEN_CHOICE_SEM(:,:,p),'ok','markerfacecolor','k','markersize',10);
    ylim([0,1])
    set(gca,'box','off','xTick',INVERSE_PRICE_MEAN,'XTickLabel', tick_labels,'yTick',[0.25,0.5,0.75 1],'FontSize',15)
    title(sprintf('%s', partners{p}),'FontSize',17)
    ylabel('p(generous)','FontSize',15); xlabel('Donation Efficiency','FontSize',15)
end

% ANOVA

Y = [GEN_CHOICE(:,:,1)]; % friend
Y = Y(:);
W = [GEN_CHOICE(:,:,2)]; % stranger
W = W(:);
Y = [Y;W];
group = combvec(1:30,1:5,1:2)';
[p,~, stats]=anovan(Y, group,'random',1,'model','interaction','varnames',{'subject','price','partner'});

% p =
% 
%    0.194816627834495 % subject
%    0.000000000000000 % price
%    0.000020013268110 % partner
%    0.000366068013740 % subject*price
%    0.000000002888552 % subejct*partner
%    0.000766272735223 % price*partner

%%

% GO TO TRIALS SELECTION STUDY 2


%% - PARAMETERS OF ALL MODELS 

n_subjects = numel(all_studies);

PARAMS_CR = nan(n_subjects,3,2);% 3=num parameters of this model; 2=num partners
PARAMS_CD = nan(n_subjects,2,2);
PARAMS_R = nan(n_subjects,2,2);

for p=1:n_partners
    for s=1:n_subjects
        PARAMS_CR(s,1,p)=AS{s}.parshatM(1,1,p);
        PARAMS_CR(s,2,p)=AS{s}.parshatM(1,2,p);
        PARAMS_CR(s,3,p)=AS{s}.parshatM(1,3,p);
        
        PARAMS_CD(s,1,p)=AS{s}.parshatM(2,1,p);
        PARAMS_CD(s,2,p)=AS{s}.parshatM(2,2,p);
        
        PARAMS_R(s,1,p)=AS{s}.parshatM(3,1,p);
        PARAMS_R(s,2,p)=AS{s}.parshatM(3,2,p);
    end
end

PARAMS{1}=PARAMS_CR;
PARAMS{2}=PARAMS_CD;
PARAMS{3}=PARAMS_R;

% save('AS_CESS','PARAMS','PARAMS_CR','PARAMS_CD','PARAMS_R','-append')



% ANOVA ON PARAMS - STUDY 2

n_subjects = numel(study_2);

Y  = [PARAMS_CR(study_2,2,1);PARAMS_CR(study_2,3,1);PARAMS_CR(study_2,2,2);PARAMS_CR(study_2,3,2)];
group = combvec(1:n_subjects,1:2,1:2)'; % subject / context/ partner
[p,~, stats]=anovan(Y, group,'random',1,'model','interaction','varnames',{'subject','status','partner'});

% without interaction
% p =
% 
%    1.0e-03 *
% 
%    0.000000000000001
%    0.000028344444225
%    0.327426737084570

% with interaction
% 
% p =
% 
%    0.000004337275146 % subject
%    0.000090220785907 % status (better or worse off)
%    0.000373061936493 % partner
%    0.000009275542432 %
%    0.006143208243300 %
%    0.073417124528679 % status*partner

% study 2:

mean_ro_f = mean (PARAMS_CR(study_2,2,1)) % = 0.3039
sem_ro_f = std (PARAMS_CR(study_2,2,1))/sqrt(n_subjects) % = 0.0663
[H,P,CI,STATS] = ttest(PARAMS_CR(study_2,2,1)) % p = 8.0591e-05 ; t = 4.583

mean_ro_s = mean(PARAMS_CR(study_2,2,2)) % = 0.126
sem_ro_s = std (PARAMS_CR(study_2,2,2))/sqrt(n_subjects) % 0.0713
[H,P,CI,STATS] = ttest(PARAMS_CR(study_2,2,2)) % p = 0.0874; t= 1.7693

mean_sig_f = mean(PARAMS_CR(study_2,3,1)) % = 0.0392
sem_sig_f = std (PARAMS_CR(study_2,3,1))/sqrt(n_subjects) % = 0.0778
[H,P,CI,STATS] = ttest(PARAMS_CR(study_2,3,1)) % p=0.6181; t=0.5040

mean_sig_s = mean(PARAMS_CR(study_2,3,2)) % = -0.0594
sem_sig_s = std (PARAMS_CR(study_2,3,2))/sqrt(n_subjects) % 0.0622
[H,P,CI,STATS] = ttest(PARAMS_CR(study_2,3,2)) % p = 0.3477; t = -0.9544


