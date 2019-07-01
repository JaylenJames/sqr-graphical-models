function [Xt, labels] = load_count_dataset(dataset, nDim)
%load_count_dataset - Loads count dataset and filters to nDim
%  With the exception of "brca", each dataset is reduced to nDim
%  by sorting the variables by the sum/mean over all instances.
%  For "brca", the variables are sorted by variance instead of 
%  sum/mean because variance is often more important in biology.

Xt = []; labels = {}; % Placeholders

if(strcmp(dataset,'brca'))
    %% If "brca", first filter variables by variance, then take log of raw counts
    % Load raw BRCA counts
    load(sprintf('../data/%s-raw.mat',dataset),'Xt','labels');
    
    % Filter based on variance
    [~,sortedI] = sort(var(log(Xt+1)),'descend');
    Xt = Xt(:,sortedI(1:nDim));
    labels = labels(sortedI(1:nDim));
    
    % Sort by sum/mean
    [Xt, labels] = transform( Xt, labels, nDim );

    % Take log(x+1) transformation
    Xt = floor(log(Xt+1));
else
    %% Otherwise simply load dataset and filter dimensions
    load(sprintf('../data/%s.mat',dataset),'Xt','labels');
    [Xt, labels] = transform( Xt, labels, nDim );
end

% Error check
if(nDim > size(Xt,2)); warning('Requested number of dimensions (nDim = %d) is larger than the loaded dataset \n which has dimension (d = %d)', nDim, size(Xt,2)); end
fprintf('\nSuccessfully loaded dataset %s with d = %d, n = %d, and sum of all counts = %d\n\n', dataset, size(Xt,2), size(Xt,1), full(sum(Xt(:))));
end

%% Sort and filter to nDim dimensions
function [ transXt, transLabels] = transform( Xt, labels, nDim)
    % Determine wordCounts
    wordCounts = full(sum(Xt));
    [~, idx] = sort(wordCounts, 'descend');

    % Sort rows based on wordCounts
    transXt = Xt(:, idx);
    transLabels = labels(idx);

    % Filter words if not -1 or negative
    if(nDim > 0 && nDim <= length(transLabels))
        transXt = transXt(:,1:nDim);
        transLabels = transLabels(1:nDim);
    end
end