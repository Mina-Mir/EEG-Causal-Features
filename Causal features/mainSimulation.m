% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [causesFoundHsic] = mainSimulation(y0, y1, r)

info.minDependencyThreshold = 0.05;
info.thresholdCI = 0.5;
info.dependencyThresholdT0T1 = 0.05;
info.methods = ["Hsic"];
info.typeDependencyRY1 = "hsic";
info.typeDependencyY1Y0 = "hsic";
% swLibraryDir = '../KCI';
% addpath(genpath(swLibraryDir));


% run level1 dependency tests
[candidatesHsicLevel1, dependencyThresholdHsic]...
    = calcCandidatesLevel1(info, y1, r);

%% find connections Y0 Y1 + conditional independence test Yo_||_R|Y1 + locate features + save results
dependencyThreshold = dependencyThresholdHsic;
candidatesLevel1 = candidatesHsicLevel1;
info.typeDependencyRY1 = "hsic";
info.typeDependencyY1Y0 = "hsic";
%% Find connections y0y1
if  strcmp(info.typeDependencyY1Y0, "hsic")
	[statHSIC, pHSIC] = calcDependencyY0Y1(y0(:,candidatesLevel1), y1(:,candidatesLevel1), info.typeDependencyY1Y0);
    strongConnectionsY0Y1 = (pHSIC < info.dependencyThresholdT0T1);
end

%% Calculate CItest on root node
pars.pairwise = true;	    % if true, the test is performed pairwise if d1>1 (standard: false)
pars.bonferroni  = false;   % if true, bonferroni correction is performed (standard: false)
pars.width = 2;             % kernel width (standard: 0, which results in an automatic -heuristic- choice)

nCandidates = length(candidatesLevel1);

pvalHsic = zeros(1, nCandidates);

parfor iCand = 1 :  nCandidates
	conditionToTheseCandidates = candidatesLevel1(iCand);         

    if strongConnectionsY0Y1(iCand) % if y^1_0->y^1_1 exists
        checkCIVar = candidatesLevel1(iCand);  % check y^1_0 _||_ R | y^1_1
        [pvalHsic(iCand), ~] = indtest_hsic(y0(:, checkCIVar), r, y1(:,conditionToTheseCandidates), pars);
    end
end

causesFoundHsic = candidatesLevel1(pvalHsic > info.thresholdCI);
