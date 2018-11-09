function [InitialVars] = structure_initial_variables (fitmodels, nfmodels, maxnumpars, nparshat_vec)

InitialVars.fitmodels = fitmodels;
InitialVars.nfmodels = nfmodels;
InitialVars.maxnumpars = maxnumpars;
InitialVars.nparshat_vec = nparshat_vec;
%InitialVars.partners = partners; 
%InitialVars.n_partners = n_partners;
%InitialVars.s_n_list = s_n_list;
%InitialVars.n_subjectcts = n_subjects;

end


