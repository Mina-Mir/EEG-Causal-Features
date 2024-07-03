% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [statHSIC, pHSIC] = calcDependencyY1R(y1, r)
% calculate dependency between M^i and R 
pars.pairwise = true;	
pars.bonferroni = false;
pars.perm = 1000;
nComponents = size(y1,2);
pHSIC = nan(1, nComponents);
statHSIC = nan(1, nComponents);

parfor iComp = 1: nComponents
    [pHSIC(iComp), statHSIC(iComp)] = indtest_hsic(y1(:,iComp), r,[],  pars);
end
end
