% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [statHSIC, pHSIC] = calcDependencyY0Y1(y0, y1, type)
% calculate dependency between P^i and M^i variables
pars.pairwise = true;	
pars.bonferroni = false;
pars.perm = 1000;
nCandidates = size(y1,2);
pHSIC = nan(1, nCandidates);
statHSIC = nan(1, nCandidates);

for iComp = 1 : nCandidates
	if strcmp(type, "hsic")
		[pHSIC(iComp), statHSIC(iComp)] = indtest_hsic(y0(:,iComp), y1(:,iComp),[],  pars);
	end
end
end
