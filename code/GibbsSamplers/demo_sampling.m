%% Load data
addpath('GibbsSamplers/');
%load('three-word-illustration.mat'); % Original fitted distribution
alpha = [0.9539, 0.4426, 0.0404]';
neg = -496;
pos = 0.198;
Theta = [0   neg neg;
         neg 0   pos;
         neg pos 0  ];
words = {'library','boundary','layer'};

%% Sample using Gibbs TPGM
n = 1000;
p = 3;
R = 20; % True maximum from is closer to 23 but clipped to 20
for maxit = [25, 100, 400, 1600]
    %% Sample with a certain maxit and show plot
    DW = GibbsTPGM(n,p,R,alpha,Theta,maxit);

    % Plot
    empiricalHist = zeros(repmat(length(0:R),1,p));
    for i = 1:n
        empiricalHist(DW(i,1)+1, DW(i,2)+1, DW(i,3)+1) = empiricalHist(DW(i,1)+1, DW(i,2)+1, DW(i,3)+1) + 1;
    end
    Zemp = (empiricalHist+eps)/sum(empiricalHist(:)); % Normalize
    ops.modelName = sprintf('TPGM Gibbs (maxit=%d, pos=%g)', maxit, pos);
    mrfs.visualizers.plotpairwisehist( Zemp, words, ops);
end

%% Plot true distribution
Z3 = zeros(repmat(length(0:R),1,p));
for i = 0:R, for j = 0:R, for k = 0:R
     x = [i,j,k]';
     Z3(i+1,j+1,k+1) = alpha'*x + x'*Theta/2*x - sum(log(factorial(x)));
end, end, end;
[~,maxI] = max(Z3(:));
[a,b,c] = ind2sub(size(Z3),maxI);
maxX = [a,b,c]'-1;
Z3 = exp(Z3 - max(Z3(:)));
Z3 = Z3/sum(Z3(:)); % Normalize
ops.modelName = sprintf('True (pos=%g, Exhaustive enumeration)',pos);
mrfs.visualizers.plotpairwisehist( Z3, words, ops );

%% Display reasoning that quadratic term overcomes base measure
maxX
linTerm = alpha'*maxX
quadTerm = maxX'*Theta/2*maxX
baseTerm = -sum(log(factorial(maxX)))

logProportion = linTerm + quadTerm + baseTerm
logProportionOfZero = 0;

ratioProb = exp(logProportion)/exp(logProportionOfZero)