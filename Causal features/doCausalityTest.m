% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function outputTable = doCausalityTest(info, y0, y1,  r, group, candidatesLevel1, dependencyThreshold)

nMethods = 1;
%% Find connections y0y1
if strcmp(info.typeDependencyY1Y0, "hsic")
	[statHSIC, pHSIC] = calcDependencyY0Y1(y0(:,candidatesLevel1), y1(:,candidatesLevel1), info.typeDependencyY1Y0);
    strongConnectionsY0Y1 = (pHSIC < info.dependencyThresholdT0T1);
end

%% Calculate CItest on root node
pars.pairwise = true;	    % if true, the test is performed pairwise if d1>1 (standard: false)
pars.bonferroni  = false;   % if true, bonferroni correction is performed (standard: false)
pars.width = 0;             % kernel width (standard: 0, which results in an automatic -heuristic- choice)

nCandidates = length(candidatesLevel1);
pval = zeros(1, nCandidates);

for iCand = 1 :  nCandidates
	conditionToTheseCandidates = candidatesLevel1(iCand);         

    if strongConnectionsY0Y1(iCand) % if y^1_0->y^1_1 exists
        checkCIVar = candidatesLevel1(iCand);  % check y^1_0 _||_ R | y^1_1
        [pval(iCand), ~] = indtest_hsic(y0(:, checkCIVar), r, y1(:,conditionToTheseCandidates), pars);
    else
        disp('Candidate rejected as there was no connection with past rootnode . Maybe falsly rejection')
    end
end

%% locate causal features

nElectrodes = length(info.electrodes);
candidatesY0Y1 = candidatesLevel1(pval >= info.thresholdCI);
bandsFound = ceil(candidatesY0Y1 / nElectrodes);
chansFound = mod(candidatesY0Y1, nElectrodes);
chansFound(chansFound == 0) = nElectrodes;


%% save results
% initialise vectors
maxElements = length(info.electrodes) * nMethods;
Method = strings(maxElements,1);
Subject = nan(maxElements,1);
Group = nan(maxElements,1);
ChannelFoundCausal = strings(maxElements,1);
Band = nan(maxElements,2);
ChannelsTested = strings(maxElements,1);
TypeRt = strings(maxElements,1);
TypeDependencyRY1 = strings(maxElements,1);
TypeDependencyY1Y0 = strings(maxElements,1);
DependencyThreshold = nan(maxElements,1);
ThresholdCI = nan(maxElements,1);
cnt = 0;

% Populate vectors with data
if ~isempty(chansFound)
    for iFeatFound = 1:length(chansFound)
        cnt = cnt + 1;
        Method(cnt) = info.methods{1};
        Subject(cnt) = info.subject;
        Group(cnt) = group;
        ChannelFoundCausal(cnt) = string(info.electrodes{chansFound(iFeatFound)});
        Band(cnt,:) = info.selectedBands(bandsFound(iFeatFound),:);
        ChannelsTested(cnt) = info.typeElectrodes;
        TypeRt(cnt) = info.typeRt;
        TypeDependencyRY1(cnt) = info.typeDependencyRY1;
        TypeDependencyY1Y0(cnt) = info.typeDependencyY1Y0;
        DependencyThreshold(cnt) = dependencyThreshold;
        ThresholdCI(cnt) = info.thresholdCI;
        
    end
end

outputTable = table(Method, Subject, Group, ChannelFoundCausal,...
           Band, ChannelsTested, TypeRt, TypeDependencyRY1, TypeDependencyY1Y0, DependencyThreshold, ...
           ThresholdCI);           
if cnt < maxElements
    outputTable(cnt + 1 : end, :) = [];
end
	
end
