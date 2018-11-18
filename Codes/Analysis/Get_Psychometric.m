function [Psychometric] = Get_Psychometric (s,initialVars,analysisVars)


load('AS_CESS.mat')

[~, nfmodels, ~, ~] = destructure_initial_variables(initialVars);
[partners,n_partners,~, ~] = destructure_analysis_variables(analysisVars);
n_models = nfmodels;

Psychometric = cell(n_partners,1);


for p = 1:n_partners
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
    
    % Categorize choices as generous or selfish. 1 or -1 indicating non-ref
    % or ref option
    selfish = yi - xi > 0;
    selfish = selfish*2-1;
    generous = yj - xj > 0;
    generous = generous*2 - 1;
    ind = find(selfish+generous); % CHECK IF I CAN ERASE THIS. 
    for k = 1:length(ind)
        if xj(ind(k)) > yj(ind(k))
            selfish(ind(k)) = 1;
        elseif yj(ind(k)) > xj(ind(k))
            generous(ind(k)) = 1;
        end
    end
    
    
    % categorize trials as both better off, both worse off, and one better
    % one worse
    context = double(rx+ry);
    n_contexts = 3; % or numel(unique(context));
    idx_context = nan(n_trials,n_contexts);
    
    for c=1:n_contexts
    switch c
        case 1
            idx_c = context == 2; % both better off
        case 2
            idx_c=context==1; % one and one
        case 3
            idx_c=context==0; % both worse off
    end
    idx_context(:,c) = idx_c;
    end
    
    % variable_1 = {'benefit of giving'};
    % variable_2 = {'cost_of_giving'};
    genx_idx = find(generous==-1); % -1 for choice x
    geny_idx = find(generous==1); % 1 for choice y
    
    
    
    % - BENEFIT OF GIVING: INVERSE PRICE OF GIVING OR DONATION EFFICIENCY

    inverse_price = nan(n_trials,1);
    inverse_price(genx_idx) = -(xj(genx_idx)-yj(genx_idx))./(xi(genx_idx) - yi(genx_idx));
    inverse_price(geny_idx) = -(yj(geny_idx)-xj(geny_idx))./(yi(geny_idx) - xi(geny_idx));
    inverse_price(inverse_price(:)==-Inf)=7; % INIFIT HERE - TRACK THIS!!!
    
    
    if s<=30 
        x_levels = [0 1.1;1.2 2.8;2.9 4.5;4.6 5.5;5.5 10]; 
    elseif s>30
         x_levels = [-5.1 -3.9; -3.1 -2; 0 1.1;1.2 2.8;2.9 4.5;4.6 5.5;5.5 10]; 
    end
    
    n_levels = size(x_levels,1);% 5 rows  for levels, and 2 columns, for beginning and ending
    
    
    level_idx_Mat = nan(n_trials,n_levels);
    
    for L=1:n_levels
        level_idx_Mat(:,L)=inverse_price>x_levels(L,1) & inverse_price<x_levels(L,2);
    end
    
    if sum(sum(level_idx_Mat)) == n_trials == 0
        display ('Attention! Error assisgning trials to proper price levels for psychometric ')
    end

    
    inverse_price_elements = cell(n_contexts,n_levels);
    inverse_price_mean = nan(n_contexts,n_levels);  
    inverse_price_sem = nan(n_contexts,n_levels);
    for c=1:n_contexts; 
        for L=1:n_levels;
            idxx = idx_context(:,c) & level_idx_Mat(:,L);
            inverse_price_elements{c,L} = inverse_price(idxx);
            inverse_price_mean(c,L) = mean(inverse_price(idxx));
            inverse_price_sem(c,L) = std(inverse_price(idxx))/sqrt(numel(inverse_price(idxx)));
        end
    end
            
    gen_choice_elements = cell(n_contexts,n_levels);
    gen_choice_mean = nan(n_contexts,n_levels);
    gen_choice_sem = nan(n_contexts,n_levels);
    for c=1:n_contexts; 
        for L=1:n_levels;
            idxx = idx_context(:,c) & level_idx_Mat(:,L);
            gen_choice_elements{c,L} = generous(idxx)==ct(idxx);
            gen_choice_mean (c,L) = mean(generous(idxx)==ct(idxx));
            gen_choice_sem (c,L) = std(generous(idxx)==ct(idxx))/sqrt(numel(generous(idxx)==ct(idxx)));
        end
    end
    
    model_prediction_elements = cell(n_contexts,n_levels,n_models);
    model_prediction_mean = nan(n_contexts,n_levels,n_models);
    model_prediction_sem = nan(n_contexts,n_levels,n_models);
    for m=1:n_models
        pars = PARAMS{m}(s,:,p);
        [p_generous]= model_predictions(xi, xj, yi, yj, m,pars,genx_idx,geny_idx);
        for c=1:n_contexts
            for L=1:n_levels
                idxx = idx_context(:,c) & level_idx_Mat(:,L);
                model_prediction_elements{c,L,m} = p_generous(idxx);
                model_prediction_mean(c,L,m) = mean(p_generous(idxx));
                model_prediction_sem(c,L,m) = std(p_generous(idxx))/sqrt(numel(p_generous(idxx)));
            end
        end
    end

    
    Psychometric{p}.Inverse_price_of_giving.by_context.inverse_price_elements = inverse_price_elements;
    Psychometric{p}.Inverse_price_of_giving.by_context.inverse_price_mean = inverse_price_mean;
    Psychometric{p}.Inverse_price_of_giving.by_context.inverse_price_sem = inverse_price_sem;
    
    Psychometric{p}.Inverse_price_of_giving.by_context.gen_choice_elements = gen_choice_elements;
    Psychometric{p}.Inverse_price_of_giving.by_context.gen_choice_mean = gen_choice_mean;
    Psychometric{p}.Inverse_price_of_giving.by_context.gen_choice_sem = gen_choice_sem;
    
    Psychometric{p}.Inverse_price_of_giving.by_context.model_prediction_elements = model_prediction_elements;
    Psychometric{p}.Inverse_price_of_giving.by_context.model_prediction_mean = model_prediction_mean;
    Psychometric{p}.Inverse_price_of_giving.by_context.model_prediction_sem = model_prediction_sem;
    
    
    % collapsed across contexts:
    collapsed.inverse_price_elements = cell(1,n_levels);
    collapsed.inverse_price_mean = nan(1,n_levels);
    collapsed.inverse_price_sem = nan(1,n_levels);
    for L=1:n_levels
        for c=1:n_contexts
            collapsed.inverse_price_elements{L} = [collapsed.inverse_price_elements{L}; inverse_price_elements{c,L}];
        end
        collapsed.inverse_price_mean(L) = mean(collapsed.inverse_price_elements{L});
        collapsed.inverse_price_sem(L) = std(collapsed.inverse_price_elements{L})/sqrt(numel(collapsed.inverse_price_elements{L}));
    end
    
    collapsed.gen_choice_elements = cell(1,n_levels);
    collapsed.gen_choice_mean = nan(1,n_levels);
    collapsed.gen_choice_sem = nan(1,n_levels);
    for L=1:n_levels
        for c=1:n_contexts
            collapsed.gen_choice_elements{L} = [collapsed.gen_choice_elements{L}; gen_choice_elements{c,L}];
        end
        collapsed.gen_choice_mean(L) = mean(collapsed.gen_choice_elements{L});
        collapsed.gen_choice_sem(L) = std(collapsed.gen_choice_elements{L})/sqrt(numel(collapsed.gen_choice_elements{L}));
    end
    
    collapsed.model_prediction_elements = cell(n_models,n_levels);
    collapsed.model_prediction_mean = nan(n_models,n_levels);
    collapsed.model_prediction_sem = nan(n_models,n_levels);
    for m=1:n_models
        for L=1:n_levels
            for c=1:n_contexts
                collapsed.model_prediction_elements{m,L} = [collapsed.model_prediction_elements{m,L}; model_prediction_elements{c,L,m}];
            end
            collapsed.model_prediction_mean(m,L) = mean(collapsed.model_prediction_elements{m,L});
            collapsed.model_prediction_sem(m,L) = std(collapsed.model_prediction_elements{m,L})/sqrt(numel(collapsed.model_prediction_elements{m,L}));
        end
    end
    
    collapsed.level_idx_Mat = level_idx_Mat;
    
    Psychometric{p}.Inverse_price_of_giving.collapsed = collapsed;

    
    
    % COST OF GIVING - PRICE OF GIVING
    
    price = nan(n_trials,1);
    price(genx_idx) = -(xi(genx_idx)-yi(genx_idx))./(xj(genx_idx) - yj(genx_idx));
    price(geny_idx) = -(yi(geny_idx)-xi(geny_idx))./(yj(geny_idx) - xj(geny_idx));    

    
    if s<=30
        x_levels = [0 0.01;0.1 0.35;0.36 0.6;0.65 1.1;1.2 3];
    elseif s>30
        x_levels = [-0.5 -0.31; -0.3 -0.01;0 0.01;0.1 0.35;0.36 0.6;0.65 1.1;1.2 3];
    end
    
    n_levels = size(x_levels,1);% 5 rows  for levels, and 2 columns, for beginning and ending    
    
    
    level_idx_Mat = nan(n_trials,n_levels);
    
    for L=1:n_levels
        level_idx_Mat(:,L)=price>=x_levels(L,1) & price<=x_levels(L,2);
    end
    
    if sum(sum(level_idx_Mat)) == n_trials == 0
        display ('Attention! Error assisgning trials to proper price levels for psychometric ')
    end
    
    price_elements = cell(n_contexts,n_levels);
    price_mean = nan(n_contexts,n_levels);  
    price_sem = nan(n_contexts,n_levels);
    for c=1:n_contexts; 
        for L=1:n_levels;
            idxx = idx_context(:,c) & level_idx_Mat(:,L);
            price_elements{c,L} = price(idxx);
            price_mean(c,L) = mean(price(idxx));
            price_sem(c,L) = std(price(idxx))/sqrt(numel(price(idxx)));
        end
    end
            
    gen_choice_elements = cell(n_contexts,n_levels);
    gen_choice_mean = nan(n_contexts,n_levels);
    gen_choice_sem = nan(n_contexts,n_levels);
    for c=1:n_contexts; 
        for L=1:n_levels;
            idxx = idx_context(:,c) & level_idx_Mat(:,L);
            gen_choice_elements{c,L} = generous(idxx)==ct(idxx);
            gen_choice_mean (c,L) = mean(generous(idxx)==ct(idxx));
            gen_choice_sem (c,L) = std(generous(idxx)==ct(idxx))/sqrt(numel(generous(idxx)==ct(idxx)));
        end
    end
    
    model_prediction_elements = cell(n_contexts,n_levels,n_models);
    model_prediction_mean = nan(n_contexts,n_levels,n_models);
    model_prediction_sem = nan(n_contexts,n_levels,n_models);
    for m=1:n_models
        pars = PARAMS{m}(s,:,p);
        [p_generous]= model_predictions(xi, xj, yi, yj, m,pars,genx_idx,geny_idx);
        for c=1:n_contexts
            for L=1:n_levels
                idxx = idx_context(:,c) & level_idx_Mat(:,L);
                model_prediction_elements{c,L,m} = p_generous(idxx);
                model_prediction_mean(c,L,m) = mean(p_generous(idxx));
                model_prediction_sem(c,L,m) = std(p_generous(idxx))/sqrt(numel(p_generous(idxx)));
            end
        end
    end

    
    Psychometric{p}.Price_of_giving.by_context.price_elements = price_elements;
    Psychometric{p}.Price_of_giving.by_context.price_mean = price_mean;
    Psychometric{p}.Price_of_giving.by_context.price_sem = price_sem;
    
    Psychometric{p}.Price_of_giving.by_context.gen_choice_elements = gen_choice_elements;
    Psychometric{p}.Price_of_giving.by_context.gen_choice_mean = gen_choice_mean;
    Psychometric{p}.Price_of_giving.by_context.gen_choice_sem = gen_choice_sem;
    
    Psychometric{p}.Price_of_giving.by_context.model_prediction_elements = model_prediction_elements;
    Psychometric{p}.Price_of_giving.by_context.model_prediction_mean = model_prediction_mean;
    Psychometric{p}.Price_of_giving.by_context.model_prediction_sem = model_prediction_sem;  

    
    % collapsed across contexts:
    clear collapsed
    collapsed.price_elements = cell(1,n_levels);
    collapsed.price_mean = nan(1,n_levels);
    collapsed.price_sem = nan(1,n_levels);
    for L=1:n_levels
        for c=1:n_contexts
            collapsed.price_elements{L} = [collapsed.price_elements{L}; price_elements{c,L}];
        end
        collapsed.price_mean(L) = mean(collapsed.price_elements{L});
        collapsed.price_sem(L) = std(collapsed.price_elements{L})/sqrt(numel(collapsed.price_elements{L}));
    end
    
    collapsed.gen_choice_elements = cell(1,n_levels);
    collapsed.gen_choice_mean = nan(1,n_levels);
    collapsed.gen_choice_sem = nan(1,n_levels);
    for L=1:n_levels
        for c=1:n_contexts
            collapsed.gen_choice_elements{L} = [collapsed.gen_choice_elements{L}; gen_choice_elements{c,L}];
        end
        collapsed.gen_choice_mean(L) = mean(collapsed.gen_choice_elements{L});
        collapsed.gen_choice_sem(L) = std(collapsed.gen_choice_elements{L})/sqrt(numel(collapsed.gen_choice_elements{L}));
    end
    
    collapsed.model_prediction_elements = cell(n_models,n_levels);
    collapsed.model_prediction_mean = nan(n_models,n_levels);
    collapsed.model_prediction_sem = nan(n_models,n_levels);
    for m=1:n_models
        for L=1:n_levels
            for c=1:n_contexts
                collapsed.model_prediction_elements{m,L} = [collapsed.model_prediction_elements{m,L}; model_prediction_elements{c,L,m}];
            end
            collapsed.model_prediction_mean(m,L) = mean(collapsed.model_prediction_elements{m,L});
            collapsed.model_prediction_sem(m,L) = std(collapsed.model_prediction_elements{m,L})/sqrt(numel(collapsed.model_prediction_elements{m,L}));
        end
    end
    Psychometric{p}.Price_of_giving.collapsed = collapsed;    
    
end

end


        
    