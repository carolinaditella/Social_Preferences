function [InitialVars] = structure_initial_variables (fitmodels, nfmodels, maxnumpars, nparshat_vec)

InitialVars.fitmodels = fitmodels;
InitialVars.nfmodels = nfmodels;
InitialVars.maxnumpars = maxnumpars;
InitialVars.nparshat_vec = nparshat_vec;

end


