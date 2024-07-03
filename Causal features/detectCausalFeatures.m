% Copyright (C)
% Max Planck Institute for Intelligent Systems,
% Atalanti A. Mastakouri atalanti.mastakouri@tuebingen.mpg.de

function [outputTableHsicHsicSub] =...
    detectCausalFeatures(info)

    outputTableHsicHsicSub = [];
    
	% Load data from mat files that contain post processed EEG data
	[y0, y1, r, group] = loadFeaturesAndResponse(info);
	if info.logRt
        r = log(r);
    end
	% run level1 dependency tests
	[candidatesHsicLevel1, dependencyThresholdHsic]	= calcCandidatesLevel1(info, y1, r);
	
	%% find connections Y0 Y1 + conditional independence test Yo_||_R|Y1 + locate features + save results

    if info.testCase4
        if ~isempty(candidatesHsicLevel1)
            dependencyThreshold = dependencyThresholdHsic;
            candidatesLevel1 = candidatesHsicLevel1;
            info.typeDependencyRY1 = "hsic";
            info.typeDependencyY1Y0 = "hsic";
            outputTableHsicHsicSub = doCausalityTest(info, y0, y1,  r, group, candidatesLevel1, dependencyThreshold);
        end
    end
	
end
