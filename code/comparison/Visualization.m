classdef Visualization
    methods (Static)
        function visMethod(cvArray, XtAll, labels, datasetLabel, filename, pMax)
            if( nargin < 4); datasetLabel = 'Unknown Dataset'; end
            if( nargin < 5); filename = []; end
            if( nargin < 6); pMax = 10; end
            XtSampleAll = cell2mat({cvArray.XtSample}');

            % Setup XtAll and XtSampleAll
            pLimit = min(pMax,size(XtAll,2));
            if(size(XtAll,2) > pMax)
                warning('Only visualizing %d dimensions even though there are %d dimensions in the dataset.',...
                    pMax, size(XtAll,2));
            end
            XtAll = full(XtAll(:,1:pLimit));
            XtSampleAll = full(XtSampleAll(:,1:pLimit));
            maxX = round(quantile(XtAll(:),0.99));

            % Clip labels
            maxLabelLength = 7;
            for li = 1:length(labels)
                if(length(labels{li}) > maxLabelLength)
                    labels{li} = [labels{li}(1:(maxLabelLength-2)), '..'];
                end
            end

            % Setup colormap
            colormap('gray');
            temp = colormap();
            colormap(flipud(temp));

            % Show marginals of true data
            maxVal = 0.4;
            subplot(1,3,1);
            poissonplotmatrix(struct('maxX',maxX,'labels',{labels},'maxVal',maxVal), XtAll);
            title(sprintf('True Marginals: %s', mapDataset(datasetLabel)));

            % Show sampled marginals
            subplot(1,3,2);
            poissonplotmatrix(struct('maxX',maxX,'labels',{labels},'maxVal',maxVal), XtSampleAll);
            title(sprintf('Sampled Marginals: %s', mapMethod(cvArray(1).methodName)));

            % Compute mean over CV folds of maxMmdPair
            meanPair = cvArray(1).pairValues;
            for cvi = 2:size(cvArray,1)
                meanPair = meanPair + cvArray(cvi).pairValues;
            end
            meanPair = meanPair./size(cvArray,1);

            % Show mmd-pair as block
            subplot(1,3,3);
            % Simple hack to display meanMaxMmdPair values
            poissonplotmatrix(struct('maxX',maxX,'labels',{labels},'meanPair',meanPair,'maxVal',maxVal), XtSampleAll);
            
            % Setup title
            if(strcmp(cvArray(1).metric,'spearman'))
                metricLabel = 'Spearman''s \rho';
            else
                metricLabel = 'MMD';
            end
            title(sprintf('Pairwise %s: %s', metricLabel, mapMethod(cvArray(1).methodName)));

            % Set position and background color
            set(gcf,'Position',[0 0 1372 437]);
            set(gcf, 'Color', 'w');
            
            % Save figure
            if(~isempty(filename))
                warning('off','all');
                export_fig(filename,'-q100');
                warning('on','all');
            end
        end
        
        function visMultHistogram(cvAll, datasetLabel, filename)
            if( nargin < 2); datasetLabel = 'Unknown Dataset'; end
            if( nargin < 3); filename = []; end
            
            methodNames = {cvAll(1,:).methodName};
            p = size(cvAll(1,1).pairValues,1);

            % Create data matrix for histogram values (each column is a method)
            histX = [];
            if(strcmp(cvAll(1).metric,'spearman'))
                % Do not include diagonal for Spearman's rho (since it's trivially 1)
                upperTri = triu(true(size(cvAll(1,1).pairValues)), 1);
            else
                % Include diagonal for MMD (since the diagonal is MMD on one dimension)
                upperTri = triu(true(size(cvAll(1,1).pairValues)), 0);
            end
            
            % Loop through each method and all CV-folds to average the pairwise values
            for mi = 1:size(cvAll,2)
                meanPairValues = cvAll(1,mi).pairValues;
                for cvi = 2:size(cvAll,1)
                    meanPairValues = meanPairValues + cvAll(cvi,mi).pairValues;
                end
                meanPairValues = meanPairValues./size(cvAll,1);
                histX = [histX, meanPairValues(upperTri)];
            end
            
            % Reorder
            [~,sortedI] = sort(mean(histX),'descend');
            histX = histX(:,sortedI);
            methodNames = methodNames(sortedI);

            if(strcmp(cvAll(1).metric,'spearman'))
                xLabelStr = 'Pairwise Spearman''s \rho Difference';
            else
                xLabelStr = 'Pairwise Maximum Mean Discrepancy';
            end
            methodNames = mapMethod(methodNames);
            multhistogram(struct('labels',{methodNames},'xlabel',xLabelStr), histX);
            title(sprintf('%s (d = %d)', mapDataset(datasetLabel), size(cvAll(1).XtSample,2) ));

            % Save figure
            if(~isempty(filename))
                if( p >= 100 )
                    set(gcf, 'Position', [0 0 474 382]);
                else
                    set(gcf, 'Position', [0 0 474 578]);
                end
                set(gcf, 'Color', 'w');

                warning('off','all');
                export_fig(filename,'-q100');
                warning('on','all');
            end
        end
    end
end

function methodCell = mapMethod(methodCell)
    isInputCell = iscell(methodCell);
    if( ~isInputCell )
        methodCell = {methodCell};
    end
    methodMap = {...
        'ind-poisson','Ind Poisson',...
        'ind-negbin','Ind Neg Bin',...
        'copula-poisson','Copula Poisson',...
        'copula-negbin','Copula Neg Bin',...
        'mixture-tune','Mixture Poiss',...
        'mixture-50','Mixture Poiss',... % Simplification for p > 100
        'ggm-tune','Gaussian GM',...
        'poisson-sqr-tune','Poisson SQR',...
        'tpgm-tune','Truncated PGM',...
        'pgm-tune','PGM',...
        'flpgm-poisson-tune','FLPGM Poisson',...
        'flpgm-negbin-tune','FLPGM Neg Bin',...
        'vine-poisson','Vine Poisson'...
        'vine-negbin','Vine Neg Bin'...
        'log-normal','Log-Normal'...
    };
    methodMap = reshape(methodMap,2,length(methodMap)/2)';
    for di = 1:length(methodCell)
        sel = strcmp(methodCell{di},methodMap(:,1));
        if(sum(sel) == 1)
            methodCell{di} = methodMap{sel,2};
        end
    end
    if( ~isInputCell )
        methodCell = methodCell{1};
    end
end

function datasetArray = mapDataset(datasetArray)
    isInputCell = iscell(datasetArray);
    if( ~isInputCell )
        datasetArray = {datasetArray};
    end
    for di = 1:length(datasetArray)
        dataset = datasetArray{di};
        switch dataset
            case '20news-with-outliers'
                dataset = '20News w/ Outliers';
            case '20news'
                dataset = '20News';
            case 'crash-severity'
                dataset = 'Crash Severity';
            case 'crime-lapd'
                dataset = 'Crime LAPD';
            case 'brca'
                dataset = 'BRCA';
            case 'classic3'
                dataset = 'Classic3';
            otherwise
                dataset = datasetArray{di}; % Default keep label
        end
        datasetArray{di} = dataset;
    end
    if( ~isInputCell )
        datasetArray = datasetArray{1};
    end
end