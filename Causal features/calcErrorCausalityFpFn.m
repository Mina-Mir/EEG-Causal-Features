% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [statResultsTable] = calcErrorCausalityFpFn
%CALCERRORCAUSALITY Summary of this function goes here
%   Detailed explanation goes here

addpath(genpath('../KCI'));
% number of features in Yo type nodes
nAll = [5 20 50 125 200];

% define probability of an edge
pAll = [0.2 0.3 0.4, 0.5];
nIter = 20;
nTrialsAll = [100 200 400 600 800 1000];
rng('default');

for iN = 1 : length(nAll)
    statResultsTable = table;
    % coefficients of linear functions f1.  f2, f3, f4
    try
        for iP = 1:length(pAll)
            p = pAll(iP);
            for iNTrials = 1 : length(nTrialsAll)
                nTrials = nTrialsAll(iNTrials);
                n = nAll(iN);
                
                falsePositivesAll =  nan(1,nIter);
                falseNegativesAll = nan(1,nIter);
                tic
                parfor iIter = 1:nIter % parfor

                    mu = 0;
                    sigma = random('uniform', 0, 1, 1);
                    gm = gmdistribution(mu,sigma);
                    parameters = random(gm, 10);
                    
                    % create noise terms for all nodes of graph
                    noiseTerms = createNodesWithNoiseTerms(n*2 + 1,nTrials); % noiseTerms = trials x nodes
                    edgeBetweenY0Y0 = binornd(1,p*ones(n));
                    edgeBetweenY0Y1 = binornd(1,p*ones(n));
                    edgeBetweenY1Y1 = binornd(1,p*ones(n));
                    edgeBetweenY0R = zeros(1,n); % allow for direct edges from P to R as well
                    edgeBetweenY1R = zeros(1,n);
                    while  sum(edgeBetweenY1R) == 0
                        edgeBetweenY1R = logical(binornd(1,p*ones(1,n)));
                        edgeBetweenY0R = logical(binornd(1,p*ones(1,n))); 
                    end
                    
                    for iNodeY0 = 1 : n
                        % remove auto-cycles
                        edgeBetweenY0Y0(iNodeY0, iNodeY0) = 0;
                        edgeBetweenY1Y1(iNodeY0, iNodeY0) = 0;
                        % remove cycles in graph (DAG)
                        edgeBetweenY0Y0(iNodeY0, [1:iNodeY0] < iNodeY0) = 0;
                        edgeBetweenY1Y1(iNodeY0, [1:iNodeY0] < iNodeY0) = 0;
                        edgeBetweenY0Y1(iNodeY0, [1:iNodeY0] < iNodeY0) = 0;
                        % enforce edge for Y^i_0 -> Y^i_1
                        edgeBetweenY0Y1(iNodeY0, iNodeY0) = 1;
                    end
                    
                    % create nodes for Y0
                    [Y0] = noiseTerms(:, 1:n);
                    % define Y0 to Y0 connections
                    Y0new = Y0;
                    for iNodeY0 = 1 : n
                        if sum(edgeBetweenY0Y0(:,iNodeY0)) >= 1
                            Y0enteringY0i = edgeBetweenY0Y0(:,iNodeY0)==1;
                            Y0new(:,iNodeY0) = Y0(:,iNodeY0) + sum(parameters(1) .* Y0(:,Y0enteringY0i) + parameters(5), 2);
                        end
                    end
                    
                    % define Y0 to Y1 connections edge always exits Y0 and enters Y1
                    [Y1] = noiseTerms(:, n+1:2*n);
                    for iNodeY1 = 1 : n
                        Y0senteringY1i = edgeBetweenY0Y1(:,iNodeY1) == 1;
                        Y1(:,iNodeY1) = Y1(:,iNodeY1) + sum(parameters(2) .* Y0new(:,Y0senteringY1i) + parameters(6), 2);
                    end
                    
                    % define Y1 to Y1 connections
                    Y1new = Y1;
                    for iNodeY1 = 1 : n
                        if sum(edgeBetweenY1Y1(:,iNodeY1)) >= 1
                            Y1enteringY1i = edgeBetweenY1Y1(:,iNodeY1)==1;
                            Y1new(:,iNodeY1) = Y1(:,iNodeY1) + sum(parameters(3) .* Y1(:,Y1enteringY1i) + parameters(7), 2);
                        end
                    end
                    
                    % define Y1 to R connections
                    [R] = noiseTerms(:, 2*n+1);
                    
                    R = R + sum(parameters(4) .* Y1new(:, edgeBetweenY1R) + parameters(8), 2) + sum(parameters(9) .* Y0new(:, edgeBetweenY0R) + parameters(10), 2);
                    
                    allConnections = [edgeBetweenY0Y0, edgeBetweenY0Y1, edgeBetweenY0R';...
                        zeros(n, n), edgeBetweenY1Y1, edgeBetweenY1R';...
                        zeros(1,2*n+1)];
                    
                    [causesFoundHsic] = mainSimulation(Y0new, Y1new, R);
                    negativesFoundHsic = setdiff([1:n], causesFoundHsic);
                    %visualizeGraph(n,allConnections);
                    [trueDirCauses, trueIndirCauses] = findCausesFromConnections(allConnections);
                    allCauses = [trueDirCauses; trueIndirCauses]';
                    trueNonCauses = setdiff([1:n], allCauses);
                    falseNegativesAll(iIter) = length(intersect(negativesFoundHsic, allCauses));
                    falsePositivesAll(iIter) = length(intersect(causesFoundHsic, trueNonCauses));
                    
                end
                
                statResultsTableN = table;
                statResultsTableN.nSamples = nTrials;
                statResultsTableN.nNodes = n;
                statResultsTableN.sparseness = p;
                statResultsTableN.nIterations = nIter;
                statResultsTableN.FPAllmean = mean(falsePositivesAll);
                statResultsTableN.FNAllmean = mean(falseNegativesAll);
                statResultsTableN.FPAllmedian = median(falsePositivesAll);
                statResultsTableN.FNAllmedian = median(falseNegativesAll);
                statResultsTableN.falsePositivesAll = falsePositivesAll;
                statResultsTableN.falseNegativesAll = falseNegativesAll;
                statResultsTableN.timeDuration = toc;
                statResultsTable = [statResultsTable; statResultsTableN];
                
                disp(strcat("sparsity: ", num2str(p), ", nSamples: ",num2str(nTrials), ", nNodes: ", num2str(n)));
                toc
            end
        end
        save(strcat('statTableCausesGraphsFPandFN_GNChangeSTDPR025_', num2str(n),'.mat'), 'statResultsTable');
    catch
        disp('error');
        save(strcat('statTableCausesGraphsFPandFN_GNChangeSTDPR025_', num2str(n),'.mat'), 'statResultsTable');
    end
end



