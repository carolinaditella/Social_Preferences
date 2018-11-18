function [mean_AEI, block_wise_AEI, temporal_AEI] = get_GARP(s_n)


% calculates AEI (Afriat's efficiency index - CCEI Critical Cost Efficiency Index) by block, 
% in chronological order. Then seprates by partner - friend or stranger


load(['/Users/carolinaditella/Dropbox/Social Preferences/OPTION PRESENTATION/SOC_PREF_DATA/',num2str(s_n),'/SOC_PREF_DATA_',num2str(s_n),'.mat']);

n_partners = 2;
    
    if s_n<=330 % subject_list_Study1(end) % define these two variables on analysisVars??
        n_b=4;
    elseif s_n>=331 % subject_list_Study2(1) % define these two variables on analysisVars??
        n_b=6;
    end
    
    AEI = nan(n_b,n_partners);
    
    for h=1:n_partners
        % CCEI with Friend:
        if h==1
            options = Data.choice_friend(:,1:4);
            ct = Data.choice_friend(:,5);
            block_info = Data.choice_friend(:,8);
        elseif h==2
            options = Data.choice_stranger(:,1:4);
            ct = Data.choice_stranger(:,5);
            block_info = Data.choice_stranger(:,8);
        end
        block_id = unique(block_info);
        %find(isnan(ct))
        options(isnan(ct),:)=[];
        ct(isnan(ct))=[];
        block_info(isnan(ct))=[];
        
        n_total_trials = size(options,1);
        
        DS = options(:,3) - options(:,1);
        DO = options(:,4) - options(:,2);
        TI = DS~=0 & DO~=0;
        TII = DS==0 & DO~=0;
        TIII = DS~=0 & DO==0;
        if sum(TI+TII+TIII)==n_total_trials %2934
            disp('categorization of trial types is alright!')
        end
        
        p = nan(size(options,1),2);
        
        % TI type trials: normal
        p(TI,1) = -DO(TI)./DS(TI);
        p(TI,2) = 1;
        
        % TII type trials: selflex trials DS==0;
        p(TII,1) = 1;
        p(TII,2) = 0;
        
        % TIII type trials: 
        p(TIII,1) = 0;
        p(TIII,2) = 1;
        
        TI(sign(p(:,1))<0)=0;
        
        % define chosen bundle:
        x = nan(size(options,1),2);
        for t=1:size(options,1)
            x (t,:)= options(t,2+ct(t):2+ct(t)+1); % 1:2 or 3:4
        end
        
        %[options ct x];
        options = options(TI,:); %options = options(1:11,:);
        ct = ct(TI); %ct = choice -1 or 1
        p = p(TI,:); %p = price
        x = x(TI,:); %x = chosen bundle
        
        for b=1:n_b
            options_temp = options(block_info(TI)==block_id(b),:);
            ct_temp = ct(block_info(TI)==block_id(b));
            p_temp = p(block_info(TI)==block_id(b),:);
            x_temp = x(block_info(TI)==block_id(b),:);
            
            V_GARP=[];
            
            O=size(p_temp,1); % this is a capital letter o, not a zero. Here it is 25, the number of trials per block
            m=eye(O);
            for i=1:O
                for j=[1:i-1, i+1:O]
                    if p_temp(i,:)*x_temp(i,:)'>=p_temp(i,:)*x_temp(j,:)' % at price i, bundle j was affordable but was not chosen. bundle i was chosen. bundle revelaed prefered to bundle j 
                        m(i,j)=1;
                    end
                end
            end
            
            mt=sign(mpower(m,(length(m))));
            for i=1:O
                for j=[1:i-1, i+1:O]
                    if mt(i,j) && p_temp(j,:)*x_temp(j,:)'>p_temp(j,:)*x_temp(i,:)' % at price j, bundle i was affordable and was not chosen, despite having been revealed prefered. Violation of rationality. 
                        V_GARP=[V_GARP;s_n i x(i,:) j x(j,:)];
                    end
                end
            end
            
            
            tic
            e=1;
            for i=1:O
                for j=[1:i-1, i+1:O]
                    while mt(i,j) && e*p_temp(j,:)*x_temp(j,:)'>p_temp(j,:)*x_temp(i,:)' % start adjusting the budget set, tiny bit by bit, until violation of rationality is erased. 
                        e=e-0.001; % 0.001
                    end
                    AEI(b,h)=e; % h is partner here. b is block. 
                end
            end
            toc
        end
    end


block_wise_AEI = AEI;
mean_AEI = mean(AEI,1); % mean by partner
temporal_AEI = block_wise_AEI(Data.block_order); % AEI in cronological order. 

end




