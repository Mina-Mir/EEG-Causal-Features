% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [candidatesHsicLevel1, dependencyThresholdHsic]...
	= calcCandidatesLevel1(info, y1, r)
  
[statHSIC, pHSIC] = calcDependencyY1R(y1, r);
% find significant dependency candidates y1_|/|_ Rt
dependencyThresholdHsic = info.minDependencyThreshold;
candidatesHsicLevel1 = find(pHSIC < dependencyThresholdHsic);

end

