% COMPILE SCORES

n_UT = 2934;
ordered_scores_MAT = nan(n_UT,3); % 3 cluster runs

for r=1:3 % 3 independent runs of the entire thing. 
    run = r;
    filefolder = strcat('Scores_',num2str(run));
    
    J = 59; % number of cluster runs to reach 2934 trials, in batches of 50. 
    all_scores = [];
    
    for j=1:J
        filename = strcat('Score_',num2str(j));
        load (fullfile('Scores', filefolder, filename))
        all_scores = [all_scores; SCORE];
    end
    
    [ordered_scores IX] = sort(all_scores, 'descend');
    ordered_scores_MAT(:,r) = ordered_scores;
    
    save(fullfile('Scores', ['SCORE_' num2str(run)]),'all_scores','ordered_scores', 'IX') % ordered by trial ID
end

save (fullfile('Scores','ordered_scores_MAT'),'ordered_scores_MAT') 

% These scores are the average log posterior of the models given the
% response. See the equations on the Supplementary Information. To this
% average log posterior a constant must be added. Since ranking doesn't change
% due to the constant, we only add it when plotting and visualizing MI
% decay.
