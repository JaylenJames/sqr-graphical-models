function [XtArray, labelsArray, cvAll, errArray] = demo_comparison_check
%% Simple test of all methods and datasets just to see if they all run
%  Note: This does not check the correctness of the methods but rather just
%  checks if they run or not.

%% Load libraries
% ******NOTE******
% The test code also depends on R being installed with the following packages:
% 1) XMRF
% 2) VineCopula
addpath('matlab-utils'); addpath('matlab-utils/export_fig'); % Auxiliary matlab utilities
addpath('comparison'); % Comparison scripts
addpath('FastMMD'); % FastMMD implementation
addpath('QUIC'); % Fast Gaussian graphical model selection
addpath('GibbsSamplers'); % Gibbs samplers for TPGM and PGM
addpath('MVPLN_Matlab'); % Multivariate log-normal Poisson mixture using Bayesian sampling
make_comparison_mex(); % Attempt to make comparison mex files if needed for FLPGM and QUIC

%% Check loading of datasets
datasetLabelArray = {'crash-severity','crime-lapd','brca', 'classic3','20news-with-outliers','20news'};
XtArray = cell(length(datasetLabelArray),1);
labelsArray = cell(length(datasetLabelArray),1);
for di = 1:length(datasetLabelArray)
    datasetLabel = datasetLabelArray{di};
    nDim = 10;
    [XtArray{di}, labelsArray{di}] = load_count_dataset(datasetLabel, nDim);
end

%% Check if methods run on very simple crash-severity dataset
methodArray = {...
    'ind-poisson','ind-negbin',...
    'ggm-tune','mixture-tune',...
    'copula-poisson','copula-negbin',...
    'vine-poisson','vine-negbin',...
    'pgm-tune','tpgm-tune',...
    'flpgm-poisson-tune','flpgm-negbin-tune',...
    'poisson-sqr-tune',...
    'log-normal'... % Computationally most expensive
    };
metric = 'spearman'; nSamples = 10; nCV = 2;
[Xt, ~] = load_count_dataset('crash-severity', 3);

% Compute evaluations for different models
cvAll = Experiment.createCvArray(nCV,length(methodArray));
successful = true;
errArray = cell(length(methodArray),1);
for i = 1:length(methodArray)
    % Compute model and evaluate
    try
        cvAll(:,i) = Experiment.test(methodArray{i}, Xt, metric, nSamples, nCV);
    catch e
        successful = false;
        errStr = sprintf('ERROR: Method %s failed\n%s (%s)\n',methodArray{i},e.message, e.identifier);
        temp = struct2cell(e.stack);
        errStr = sprintf('%sFile: %s, function: %s, line: %d\n\n', errStr,temp{:});
        fprintf('%s',errStr);
        errArray{i} = errStr;
    end
end
if(successful)
    fprintf('\n\nAll methods sucessfully ran! :-)\n');
else
    fprintf('Some methods failed\n');
    fprintf('%s',errArray{:});
end
