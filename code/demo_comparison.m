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

%% Setup
% Load dataset, possible datasets:
%   'crash-severity','crime-lapd','brca'
%   'classic3','20news-with-outliers','20news'
datasetLabel = 'classic3';
nDim = 10;
[Xt, labels] = load_count_dataset(datasetLabel, nDim);

% Setup parameters for experiment
metric = 'mmd'; % Either 'mmd' or 'spearman'
nSamples = 1000;
nCV = 3;

%% Compute evaluations for different models
% Possible methodNames:
%   'ind-poisson','ind-negbin','ggm-tune',
%   'mixture-tune','log-normal',
%   'copula-poisson','copula-negbin',
%   'vine-poisson','vine-negbin',   
%   'pgm-tune','tpgm-tune',
%   'flpgm-poisson-tune','flpgm-negbin-tune',
%   'poisson-sqr-tune'
batch = {'ind-poisson','ind-negbin'};
cvAll = Experiment.createCvArray(nCV,length(batch));
for i = 1:length(batch)
    % Compute model and evaluate
    cvAll(:,i) = Experiment.test(batch{i}, Xt, metric, nSamples, nCV);
    % Note: For methods that require hyperparameter tuning (i.e. with suffix "-tune")
    %  the default hyperparameters can be overridden by providing a parameter vector.
    %  For the mixture-tune method, the hyperparameter is the value of k.  For all other
    %  tune methods, it is the regularization value (usually denoted as lambda).
    %  For example, to attempt lambda = [1,0.1,0.01]
    %   cvAll(:,i) = Experiment.test(batch{i}, Xt, metric, nSamples, nCV, [1, 0.1, 0.01]);
    %  To merely train with a specific hyperparameter value just provide a scalar value
    %  though this is not recommended unless a hyperparameter has already been selected previously:
    %   cvAll(:,i) = Experiment.test(batch{i}, Xt, metric, nSamples, nCV, 0.01);
end

%% Visualize each method
for i = 1:length(batch)
    % Visualize method
    figure;
    filename = []; % Add filename to export figure to file
    Visualization.visMethod(cvAll(:,i), Xt, labels, datasetLabel, filename)
end

%% Visually compare models
figure;
filename = []; % Add filename to export figure to file
Visualization.visMultHistogram(cvAll, datasetLabel, filename);
