function [model, poissMean, pVec, clusterVec] = fitpoissmix( Xt, k, alpha, maxIter )
if(nargin < 3); alpha = 1/size(Xt,1); end % 1/n
if(nargin < 4); maxIter = 100; end

% K-means initialization (10 different initializations)
nReps = 10;
bestSumSquared = Inf;
for i = 1:nReps
    rng(i);
    curClusterVec = mrfs.utils.litekmeans(Xt', k)';
    sumSquared = 0;
    for j = 1:k
        Xj = Xt(curClusterVec == j,:);
        meanJ = mean(Xj);
        sqDist = dot(meanJ,meanJ) - 2*full(Xj*meanJ') + sum(Xj.*Xj,2); % Squared distance
        sumSquared = sumSquared + sum(sqDist);
    end
    if(sumSquared < bestSumSquared)
        bestSumSquared = sumSquared;
        clusterVec = curClusterVec;
    end
end

%% EM algorithm
[n,p] = size(Xt);
poissMean = NaN(p,k);
logProb = NaN(n,k);
gammalnXt = gammaln(Xt+1);
logpoisspdf = @(mean, Xt, gammalnXt) bsxfun(@plus, bsxfun(@times, Xt, log(mean)) - gammalnXt, -mean);
for i = 1:maxIter
    %% Expectation step
    parfor j = 1:k
        Xj = Xt(clusterVec == j,:);
        if(~isempty(Xj))
            poissMean(:,j) = mean(Xj+alpha,1);
        else
            poissMean(:,j) = Xt(randi(n,1,1),:) + alpha; % Randomly assign at least one point if empty cluster
        end
    end
    
    %% Maximization step
    parfor j = 1:k
        logProb(:,j) = sum( logpoisspdf(poissMean(:,j)', Xt, gammalnXt), 2); % Independence assumed
    end
    oldClusterVec = clusterVec;
    [~,clusterVec] = max(logProb,[],2);
    
    if(all(clusterVec == oldClusterVec))
        break;
    end
end

poissMean = poissMean';
pVec = zeros(k,1);
for j = 1:k
    pVec(j) = sum(clusterVec == j);
end
pVec = pVec./sum(pVec);

model.poissMean = poissMean;
model.pVec = pVec;