    function [Psychometric] = Get_Psychometric_status (s,p,analysisVars)


load('AS_CESS.mat')

[partners,~,~,~] = destructure_analysis_variables(analysisVars);

n_models = 3;

%Psychometric


partner = partners{p};
Data = AS{s}.D{p}; % partners clean data
n_trials = size(Data,1);

xi = Data(:,1);
xj = Data(:,2);
yi = Data(:,3);
yj = Data(:,4);
ct = Data(:,5);

rx = (xi-xj)>=0;
zx = (xj-xi)>0;
ry = (yi-yj)>=0;
zy = (yj-yi)>0;

% variable_1 = {'efficiency of donation'};
% variable_2 = {'cost of donation'};

rel_status_A = yi-yj; % Non-Ref = 1
rel_status_B = xi-xj; % Ref = -1
delta_status = rel_status_A - rel_status_B;

delta_status_levels = [-80 -50; -50 -40; -35 -25; -25 -6; -5 -2; 2 5; 5 10; 11 17;40 50; 55 65; 70 80; 89 120];

n_levels = size(delta_status_levels,1);

% % % % genx_idx = find(generous==-1); % -1 for choice x
% % % % geny_idx = find(generous==1); % 1 for choice y

level_idx_Mat = nan(n_trials,n_levels);
%  0 & 1 matrix, with one 1 and all the rest 0s on every row. Each row
%  reresents a trial. The one on the row indicates to what level (columns)
%  the trial belongs. 

for L=1:n_levels
    level_idx_Mat(:,L)=delta_status>delta_status_levels(L,1) & delta_status<delta_status_levels(L,2);
end

A_choice = ct == 1;

d_status_elements = cell(1,n_levels);
d_status_mean = nan(1,n_levels);

A_choice_elements = cell(1,n_levels);
A_choice_mean = nan(1,n_levels);
% A_choice_sem = nan(1,n_levels);

model_prediction_elements = cell(n_models,n_levels);
model_prediction_mean = nan(n_models,n_levels);

for L=1:n_levels;
    idxx = level_idx_Mat(:,L);
    d_status_elements{1,L} = delta_status(find(idxx));
    d_status_mean(1,L) = mean(d_status_elements{1,L});
    
    A_choice_elements{1,L} = A_choice(find(idxx));
    A_choice_mean(1,L) = mean(A_choice_elements{1,L});
    %A_choice_sem(1,L) = mean(A_choice_elements{1,L});
    
    for m=1:n_models
        pars = PARAMS{m}(s,:,p);
        model_prediction_elements{m,L} = model_predictions_status(xi, xj, yi, yj, m,pars,idxx);
        model_prediction_mean(m,L) = mean(model_prediction_elements{m,L});
    end
end



Psychometric.d_status_elements = d_status_elements;
Psychometric.d_status_mean = d_status_mean;
Psychometric.A_choice_elements = A_choice_elements;
Psychometric.A_choice_mean = A_choice_mean;
Psychometric.model_prediction_elements = model_prediction_elements;
Psychometric.model_prediction_mean = model_prediction_mean;




end



