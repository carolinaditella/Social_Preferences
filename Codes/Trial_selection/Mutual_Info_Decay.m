% MUTUAL INFORMATION DECAY

cd('/Users/carolinaditella/Dropbox/Social Preferences/Trial Selection - Study 1')

load ('SCORES.mat')

figure
plot(1:2934, ordered_scores+log(3),'k','LineWidth', 2); hold on
plot([100 100],[0 0.25],'Color',[160 160 160]/255); %, 'LineStyle','--'
title('Mutual information decay')
xlabel('Ranked trials')
ylabel('Mutual information')
ylim([0 0.25])
set(gca,'fontsize', 15)


mean(ordered_scores(1:100)) + log(3)


% log(3) is a constant on the MI equation (See Supplementary Information)
% and represents the assumption of equal priors for each of the 3 models
% considered. 
