classdef Experiment
    
    properties (Constant)
        N_LAM = 10;
    end
    
    methods (Static)
        function cvArray = test(methodName, XtAll, metric, nSamples, nCV, tuneParamVecOverride)
            if(nargin < 3); metric = 'mmd'; end
            if(nargin < 4); nSamples = 1000; end
            if(nargin < 5); nCV = 3; end
            if(nargin < 6); tuneParamVecOverride = []; end
            
            %% Handle the difference between pairwise MMD and pairwise Spearman's rho 
            if(strcmp(metric,'spearman')); isCorr = true;
            elseif(strcmp(metric,'mmd')); isCorr = false;
            else error('metric parameter must be either ''mmd'' or ''spearman''');
            end
            
            if(isCorr)
                evalFunc = @paircorr;
            else
                evalFunc = @pairmmd;
            end
            
            %% Setup
            rng default; % For reproducibility
            [n,p] = size(XtAll); % Extract dimensions from dataset
            
            % Attempt to start parallel workers
            try
                [~,nStr] = system('grep --count "^processor" /proc/cpuinfo');
                nWorkers = str2double(nStr);
                if(isnan(nWorkers)); nWorkers = 4; end
            catch
                nWorkers = 4;
            end
            initParallel(nWorkers);
            
            % Split into CV train and test sets
            cvArray = Experiment.createCvArray(nCV,1);
            rndSeed = 1; % For reproducibility
            [trainIdxArray, testIdxArray] = mrfs.utils.cvsplit( n, nCV, rndSeed );
            sigmaVec = 10.^(-2:0.2:2)';
            nBasisPair = 2^6;% Number of basis for MMD approximation
            
            %% Run cross validation splits (most models have internal parallel programming)
            %parfor cvi = 1:nCV
            for cvi = 1:nCV
                %% Setup train/test split for this CV split
                XtTrain = XtAll(trainIdxArray{cvi},:);
                XtTest = XtAll(testIdxArray{cvi},:);
                
                lamMaxXX = full((1/size(XtTrain,1))*max(max(abs(triu(XtTrain'*XtTrain,1)))));
                [nTrain, pTrain] = size(XtTrain);
                
                %% Setup and modify parameters as needed 
                modMethodName = methodName;
                copulaType = 'Gaussian';
                modLParam = 1.5; % Parameter for FLPGM
                lamScaleVec = logspace(0,-4, Experiment.N_LAM);
                
                % Switch useNegBin for FLPGM
                useNegBin = false;
                if(strfind(modMethodName, 'flpgm'))
                    switch modMethodName(length('flpgm- '):end)
                        case 'negbin-tune'
                            useNegBin = true;
                            fprintf('Using negative binomial on length for FLPGM\n');
                            modMethodName = 'flpgm-tune';
                        case 'poisson-tune'
                            useNegBin = false;
                            fprintf('Using Poisson on length for FLPGM\n');
                            modMethodName = 'flpgm-tune';
                        otherwise
                            error('Unknown flpgm model "%s"',methodName);
                    end
                end
                
                %% Fit model and get samples
                tuneTime = []; % Default unless reset later
                fprintf('\n\n<< Starting CV = %d for model %s >>\n', cvi, methodName );
                switch modMethodName
                    case 'log-normal'
                        %% Fit model
                        ts = tic;
                        maxIter = 1000; burnin = round(0.4*maxIter);
                        [LN_mean, LN_sigma, result_beta, result_sigma, result_LL, maxLLInd] ...
                            = MVPLN_wrapper( XtTrain, maxIter, burnin );
                        trainTime = toc(ts);
                        cvArray(cvi).model.LN_mean = LN_mean;
                        cvArray(cvi).model.LN_sigma = LN_sigma;
                        cvArray(cvi).model.maxIter = maxIter;
                        cvArray(cvi).model.burnin = burnin;
                        
                        %% Visualize parameters
                        betaMat = squeeze(result_beta(:,1,:));
                        sigmaMat = result_sigma(:,1:9);
                        subplot(3,1,1);
                        plot(result_LL);
                        hold on; plot(maxLLInd*ones(2,1), ylim(),'-k'); plot(burnin*ones(2,1), ylim(),'--k'); hold off;
                        
                        subplot(3,1,2);
                        plot(betaMat);
                        xl = xlim();
                        hold on; plot([burnin,xl(2)],ones(2,1)*LN_mean, '--'); plot(maxLLInd*ones(2,1), ylim(),'-k'); plot(burnin*ones(2,1), ylim(),'--k'); hold off;
                        
                        subplot(3,1,3);
                        plot(sigmaMat);
                        xl = xlim();
                        hold on; plot([burnin,xl(2)],ones(2,1)*LN_sigma(:)', '--'); plot(maxLLInd*ones(2,1), ylim(),'-k'); plot(burnin*ones(2,1), ylim(),'--k'); hold off;
                        
                        %% Sample
                        % Sample from normal then convert to log normal
                        ts = tic;
                        ZtNormal = mvnrnd(LN_mean, LN_sigma, nSamples);
                        ZtLogNormal = exp(ZtNormal);
                        % Sample from Poisson
                        XtSample = poissrnd(ZtLogNormal);
                        sampleTime = toc(ts);
                    
                    case 'ind-poisson'
                        %% Super simple independent Poisson estimation
                        ts = tic;
                        poissMean = full(mean(XtTrain));
                        trainTime = toc(ts);
                        cvArray(cvi).model.poissMean = poissMean;
                        
                        ts = tic;
                        XtSample = poissrnd(repmat(poissMean,size(XtTrain,1),1));
                        sampleTime = toc(ts);
                        
                    case 'ind-negbin'
                        %% Fit independent negbin
                        ts = tic;
                        XtSample = NaN(nSamples,p);
                        cvArray(cvi).model.dispersion = NaN(p,1);
                        cvArray(cvi).model.poissMean = NaN(p,1);
                        cvArray(cvi).model.negbinR = NaN(p,1);
                        cvArray(cvi).model.negbinP = NaN(p,1);
                        for s = 1:p
                            xs = full(XtTrain(:,s));
                            cvArray(cvi).model.dispersion(s) = var(xs)/mean(xs);
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                params = nbinfit(xs);
                                cvArray(cvi).model.negbinR(s) = params(1);
                                cvArray(cvi).model.negbinP(s) = params(2);
                            else
                                cvArray(cvi).model.poissMean(s) = mean(xs);
                            end
                        end
                        trainTime = toc(ts);
                        
                        ts = tic;
                        for s = 1:p
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                XtSample(:,s) = nbinrnd(cvArray(cvi).model.negbinR(s), cvArray(cvi).model.negbinP(s), nSamples, 1);
                            else
                                XtSample(:,s) = poissrnd(cvArray(cvi).model.poissMean(s), nSamples,1);
                            end
                        end
                        sampleTime = toc(ts);
                    case 'copula-poisson'
                        %% Use copula-based model with Poisson marginals
                        ts = tic;
                        poissMean = full(mean(XtTrain)); % MLE estimate
                        % DT transform (take mean of DT)
                        U = (poisscdf(XtTrain, repmat(poissMean,size(XtTrain,1),1)) + poisscdf(XtTrain - 1, repmat(poissMean,size(XtTrain,1),1)))/2;
                        U(U==1) = 1-eps; % Need to make strictly < 1 because of rounding error

                        % Fit copula
                        if(strcmp(copulaType,'t'))
                            [rhohat, nuhat] = copulafit(copulaType, U);
                        else
                            rhohat = copulafit(copulaType, U);
                            nuhat = NaN;
                        end
                        trainTime = toc(ts);

                        cvArray(cvi).model.poissMean = poissMean;
                        cvArray(cvi).model.rhohat = rhohat;
                        cvArray(cvi).model.nuhat = nuhat;

                        % Get samples by getting uniform samples and then using inverse of kernel density estimate
                        ts = tic;
                        if(strcmp(copulaType, 't'))
                            USample = copularnd(copulaType, rhohat, nuhat, nSamples);
                        else
                            USample = copularnd(copulaType, rhohat, nSamples);
                        end
                        XtSample = poissinv(USample, repmat(poissMean,nSamples,1));
                        sampleTime = toc(ts);

                    case 'vine-poisson'
                        %% Use copula-based model with Poisson marginals
                        ts = tic;
                        poissMean = full(mean(XtTrain)); % MLE estimate
                        % DT transform (take mean of DT)
                        U = (poisscdf(XtTrain, repmat(poissMean,size(XtTrain,1),1)) + poisscdf(XtTrain - 1, repmat(poissMean,size(XtTrain,1),1)))/2;
                        U(U==1) = 1-eps; % Need to make strictly < 1 because of rounding error

                        %% Fit copula and sample
                        vineType = 0; indTest = 0; pLevel = 0.05; nThreads = 2; familySet = 1:6; % (i.e. try all families)
                        [vineMat, familyMat, parMat, par2Mat, USample] = ...
                            vinecopulawrapper( U, nSamples, vineType, indTest, pLevel, nThreads, familySet);
                        trainTime = toc(ts);

                        cvArray(cvi).model.vineMat = vineMat;
                        cvArray(cvi).model.familyMat = familyMat;
                        cvArray(cvi).model.parMat = parMat;
                        cvArray(cvi).model.par2Mat = par2Mat;
                        
                        %% Create actual samples using USample from vinecopulawrapper
                        ts = tic;
                        XtSample = poissinv(USample, repmat(poissMean,nSamples,1));
                        sampleTime = toc(ts);
                        
                    case 'copula-negbin'
                        %% Use copula-based model with negative binomial marginals
                        ts = tic;
                        nbinParams = NaN(p,2);
                        U = NaN(size(XtTrain));
                        cvArray(cvi).model.dispersion = NaN(p,1);
                        cvArray(cvi).model.poissMean = NaN(p,1);
                        cvArray(cvi).model.negbinR = NaN(p,1);
                        cvArray(cvi).model.negbinP = NaN(p,1);
                        for s = 1:p
                            xs = full(XtTrain(:,s));
                            cvArray(cvi).model.dispersion(s) = var(xs)/mean(xs);
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                nbinParams(s,:) = nbinfit(xs);
                                cvArray(cvi).model.negbinR(s) = nbinParams(s,1);
                                cvArray(cvi).model.negbinP(s) = nbinParams(s,2);
                                %DT Transform
                                U(:,s) = (nbincdf(xs, nbinParams(s,1), nbinParams(s,2))+nbincdf(xs-1, nbinParams(s,1), nbinParams(s,2)))/2;
                            else
                                cvArray(cvi).model.poissMean(s) = mean(xs);
                                %DT Transform
                                U(:,s) = (poisscdf(xs, cvArray(cvi).model.poissMean(s))+poisscdf(xs-1, cvArray(cvi).model.poissMean(s)))/2;
                            end
                        end
                        U(U==1) = 1-eps; % Need to make strictly < 1 because of rounding error

                        % Fit copula
                        if(strcmp(copulaType,'t'))
                            [rhohat, nuhat] = copulafit(copulaType, U);
                        else
                            rhohat = copulafit(copulaType, U);
                            nuhat = NaN;
                        end
                        trainTime = toc(ts);

                        cvArray(cvi).model.rhohat = rhohat;
                        cvArray(cvi).model.nuhat = nuhat;

                        % Get samples by getting uniform samples and then using inverse of kernel density estimate
                        ts = tic;
                        if(strcmp(copulaType, 't'))
                            USample = copularnd(copulaType, rhohat, nuhat, nSamples);
                        else
                            USample = copularnd(copulaType, rhohat, nSamples);
                        end
                        for s = 1:p
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                XtSample(:,s) = nbininv(USample(:,s), repmat(cvArray(cvi).model.negbinR(s),nSamples,1), repmat(cvArray(cvi).model.negbinP(s),nSamples,1));
                            else
                                XtSample(:,s) = poissinv(USample(:,s), repmat(cvArray(cvi).model.poissMean(s),nSamples,1));
                            end
                        end
                        sampleTime = toc(ts);

                    case 'vine-negbin'
                        %% Use copula-based model with negative binomial marginals
                        ts = tic;
                        nbinParams = NaN(p,2);
                        U = NaN(size(XtTrain));
                        cvArray(cvi).model.dispersion = NaN(p,1);
                        cvArray(cvi).model.poissMean = NaN(p,1);
                        cvArray(cvi).model.negbinR = NaN(p,1);
                        cvArray(cvi).model.negbinP = NaN(p,1);
                        for s = 1:p
                            xs = full(XtTrain(:,s));
                            cvArray(cvi).model.dispersion(s) = var(xs)/mean(xs);
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                nbinParams(s,:) = nbinfit(xs);
                                cvArray(cvi).model.negbinR(s) = nbinParams(s,1);
                                cvArray(cvi).model.negbinP(s) = nbinParams(s,2);
                                %DT Transform
                                U(:,s) = (nbincdf(xs, nbinParams(s,1), nbinParams(s,2))+nbincdf(xs-1, nbinParams(s,1), nbinParams(s,2)))/2;
                            else
                                cvArray(cvi).model.poissMean(s) = mean(xs);
                                %DT Transform
                                U(:,s) = (poisscdf(xs, cvArray(cvi).model.poissMean(s))+poisscdf(xs-1, cvArray(cvi).model.poissMean(s)))/2;
                            end
                        end
                        U(U==1) = 1-eps; % Need to make strictly < 1 because of rounding error

                        % Fit copula
                        vineType = 0; indTest = 0; pLevel = 0.05; nThreads = 2; familySet = 1:6; % (i.e. try all families)
                        [vineMat, familyMat, parMat, par2Mat, USample] = ...
                            vinecopulawrapper( U, nSamples, vineType, indTest, pLevel, nThreads, familySet);
                        trainTime = toc(ts);
                        
                        cvArray(cvi).model.vineMat = vineMat;
                        cvArray(cvi).model.familyMat = familyMat;
                        cvArray(cvi).model.parMat = parMat;
                        cvArray(cvi).model.par2Mat = par2Mat;

                        % Get samples by getting uniform samples and then using inverse of kernel density estimate
                        ts = tic;
                        for s = 1:p
                            if(cvArray(cvi).model.dispersion(s) > 1)
                                XtSample(:,s) = nbininv(USample(:,s), repmat(cvArray(cvi).model.negbinR(s),nSamples,1), repmat(cvArray(cvi).model.negbinP(s),nSamples,1));
                            else
                                XtSample(:,s) = poissinv(USample(:,s), repmat(cvArray(cvi).model.poissMean(s),nSamples,1));
                            end
                        end
                        sampleTime = toc(ts);
                        
                    case 'mixture-tune'
                        %% Setup functions
                        paramVec = 10:10:100;
                        trainFunc = @( Xt, param, model) fitpoissmix( XtTrain, param );
                        sampleFunc = @( model, nSamples ) wrapperSampleMixture( model, nSamples );
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model.poissMean = model.poissMean;
                        cvArray(cvi).model.pVec = model.pVec;
                        cvArray(cvi).model.k = fParam;
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;
                    case 'ggm-tune'
                        %% Setup functions
                        paramVec = lamScaleVec * lamMaxXX;
                        trainFunc = @( Xt, lam, model ) struct('meanX', mean(Xt),...
                            'precMat', QUIC('default', full(cov(Xt)), lam, 1e-6, 1, 100) );
                        sampleFunc = @( model, nSamples ) mvnrnd(model.meanX, inv(model.precMat), nSamples);
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model.precMat = model.precMat;
                        cvArray(cvi).model.meanX = model.meanX;
                        cvArray(cvi).model.lam = fParam;
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;
                        
                    case 'pgm-tune'
                        %% Setup functions
                        paramVec = lamScaleVec * lamMaxXX;
                        %nThreads = 20;
                        %trainFunc = @( Xt, lam, model) xmrfwrapper( Xt, 'pgm', lam, nThreads );
                        trainFunc = @( Xt, lam, model) wrapperPgm( Xt, lam );
                        
                        nGibbs = 5000;
                        sampleFunc = @( model, nSamples ) ...
                            GibbsPGM(nSamples, length(model.thetaNode), model.thetaNode, model.thetaEdge, nGibbs);
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model.thetaNode = model.thetaNode;
                        cvArray(cvi).model.thetaEdge = model.thetaEdge;
                        cvArray(cvi).model.lam = fParam;
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;
                        
                    case 'tpgm-tune'
                        %% Setup functions
                        paramVec = lamScaleVec * lamMaxXX;
                        R = full(round(quantile(XtTrain(XtTrain>0),0.99)));
                        nThreads = 20;
                        trainFunc = @( Xt, lam, model) xmrfwrapper( Xt, 'tpgm', lam, nThreads, R );
                        
                        nGibbs = 5000;
                        sampleFunc = @( model, nSamples ) ...
                            GibbsTPGM(nSamples, length(model.thetaNode), R, model.thetaNode, model.thetaEdge, nGibbs);
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model.thetaNode = model.thetaNode;
                        cvArray(cvi).model.thetaEdge = model.thetaEdge;
                        cvArray(cvi).model.lam = fParam;
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;

                    case 'flpgm-tune'
                        %% Setup functions
                        paramVec = lamScaleVec * lamMaxXX;
                        
                        trainFunc = @(Xt, lam, model) wrapperFlpgm( Xt, lam, modLParam );
                        nAnneal = 100; nGibbs = 1;
                        sampleFunc = @(model, nSamples) wrapperSamplerFlgpm( model, nSamples, nAnneal, nGibbs, useNegBin );
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model = model;
                        cvArray(cvi).model.lam = fParam;
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;
                        
                    case 'poisson-sqr-tune'
                        %% Setup functions
                        paramVec = lamScaleVec * sqrt(lamMaxXX);
                        %iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
                        %trainFunc = @(Xt, param, grm) iif( ...
                        %    nargin < 3 || isempty(grm),  @() mrfs.grm.fitgrm( Xt, 2, struct('lam', param) ), ...
                        %    true,  @() mrfs.grm.fitgrm( Xt, 2, struct('lam', param, 'grm', grm ) )...
                        %    );
                        trainFunc = @(Xt, param, grm) mrfs.grm.fitgrm( Xt, 2, struct('lam', param ) );
                        
                        nGibbs = 5000; nInner = 2; nBatches = 20;
                        sampleFunc = @( model, nSamples ) wrapperSampleSQR( model.getPsiExt(), nSamples, nGibbs, nInner, nBatches );
                        
                        %% Call tuning function
                        [model, XtSample, fParam, timing, ~] = ...
                            tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride);
                        
                        %% Save parameters
                        cvArray(cvi).model.PsiExt = model.getPsiExt();
                        cvArray(cvi).model.lam = fParam;
                        cvArray(cvi).model.lamMax = sqrt(lamMaxXX);
                        trainTime = timing.train; tuneTime = timing.tune;
                        sampleTime = timing.sample;

                    otherwise
                        error('Model %s not implemented (did you forget to add the suffix ''-tune''\n to the methods that require tuning?\n', methodName);
                end
                
                %% Compute test statistics
                pairValues = evalFunc(XtTest, XtSample, sigmaVec, nBasisPair);

                % Save some values
                cvArray(cvi).methodName = methodName;
                if(nnz(XtSample)/numel(XtSample) < 0.1)
                    cvArray(cvi).XtSample = sparse(XtSample);
                else
                    cvArray(cvi).XtSample = XtSample;
                end
                if(isCorr)
                    cvArray(cvi).evalParams = struct();
                else
                    cvArray(cvi).evalParams = struct('sigmaVec',sigmaVec,'nBasisPair',nBasisPair);
                end
                cvArray(cvi).tuneTime = tuneTime;
                cvArray(cvi).trainTime = trainTime;
                cvArray(cvi).sampleTime = sampleTime;
                cvArray(cvi).pairValues = pairValues;
                cvArray(cvi).metric = metric;
                
                % Print out some results
                if(strcmp(metric,'spearman')); upperTri = triu(true(size(pairValues)), 1);
                else upperTri = triu(true(size(pairValues)), 0); end
                fiveNum = num2cell(quantile(pairValues(upperTri),[0,0.25,0.5,0.75,1]));
                meanPair = mean(pairValues(upperTri));
                fprintf('  Model = %s\n', methodName);
                fprintf('  Metric = %s\n', metric);
                fprintf('  Mean Pair Value = %.4g\n', meanPair);
                fprintf('  Min Q1 Med Q2 Max of Pair Values = \n      [%.4g, %.4g, %.4g, %.4g, %.4g]\n', fiveNum{:});
                fprintf('  Total hyperparameter tuning time = %g s\n',trainTime);
                fprintf('  Train time = %g s\n',trainTime);
                fprintf('  Sample time = %g s\n', sampleTime);
                fprintf('<< Finished CV = %d >>\n', cvi);
            end
        end
        
        function cvArray = createCvArray(M,N)
            cvArray = struct('methodName', repmat({[]},M,N), 'model', [], 'XtSample',[],...
                'evalParams',[], 'tuneTime', [], 'trainTime',[],'sampleTime',[],'pairValues',[],'metric',[]);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pairwise evaluation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxMmdPair = pairmmd(XtTest, XtSample, sigmaVec, nBasis)
    % Compute mmd pairs in parallel
    p = size(XtTest,2);
    % Diagonals and pairwise combinations
    C = [[(1:p);(1:p)]'; nchoosek(1:p,2)];
    maxMmdPairVec = NaN(size(C,1),1);
    nC = size(C,1);
    parfor ii = 1:nC
        tts = tic;
        temp = C(ii,:);
        s = temp(1);
        t = temp(2);
        if(s == t)
            ind = s;
        else
            ind = [s,t];
        end
        [~, mmdPair] = MMDFourierFeature(XtTest(:,ind), XtSample(:,ind), sigmaVec, nBasis);
        maxMmdPairVec(ii) = max(mmdPair);
        %fprintf('  Finished pair = [%d, %d] with value = %g in %g s\n', s, t, maxMmdPairVec(ii), toc(tts));
    end
    
    % Make symmetric matrix from the values of parallel loop
    maxMmdPair = NaN(p,p);
    ind = sub2ind(size(maxMmdPair),C(:,1),C(:,2));
    maxMmdPair(ind) = maxMmdPairVec;
    ind = sub2ind(size(maxMmdPair),C(:,2),C(:,1));
    maxMmdPair(ind) = maxMmdPairVec;
end

function corrDiffPair = paircorr(XtTest, XtSample, varargin)
    corrArgs = {'type','Spearman'};
    rhoTest = corr(XtTest, corrArgs{:});
    rhoSample = corr(XtSample, corrArgs{:});
    corrDiffPair = abs(rhoTest-rhoSample);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning function for hyperparameters of model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fModel, fXtSample, fParam, timing, meanPairTune] = tuneFunc(trainFunc, paramVec, sampleFunc, XtTrain, nSamples, nCV, sigmaVec, nBasisPair, evalFunc, tuneParamVecOverride)
    if(~isempty(tuneParamVecOverride))
        paramVec = tuneParamVecOverride;
    end
    %% Tune split
    timing = [];
    ts = tic;
    rndSeed = 1;
    [XtTrainTune, XtTestTune] = mrfs.utils.traintestsplit( XtTrain, 1/nCV, rndSeed );

    %% Loop through parameter values
    fprintf('  << Starting hyperparameter tuning >>\n');
    if(length(paramVec) == 1)
        fParam = paramVec;
        bestTuneModel = [];
        meanPairTune = NaN;
        timing.tune = [];
        fprintf('    Only one parameter so not tuning\n');
    else
        meanPairTune = NaN(size(paramVec));
        medianPairTune = NaN(size(paramVec));
        nSamplesTune = max(100, nSamples/10);
        model = [];
        timing.tune.model = zeros(size(paramVec));
        timing.tune.sample = zeros(size(paramVec));
        for pi = 1:length(paramVec)
            %% Train
            tsModel = tic;
            param = paramVec(pi);
            model = trainFunc(XtTrainTune, param, model);
            timing.tune.model(pi) = toc(tsModel);

            %% Sample
            tsSample = tic;
            XtSampleTune = sampleFunc( model, nSamplesTune );
            timing.tune.sample(pi) = toc(tsSample);

            %% Compute MMD
            %[~, mmdTune] = MMDFourierFeature(XtTestTune, XtSampleTune, sigmaVec, nBasis);

            pairValues = evalFunc(XtTestTune, XtSampleTune, sigmaVec, nBasisPair);

            upperTri = triu(true(size(pairValues)),0);
            meanPairTune(pi) = mean(pairValues(upperTri));
            medianPairTune(pi) = median(pairValues(upperTri));
            fprintf('    Finished tune paramIdx = %d, param = %g, meanPair = %g, medianPair = %g in %g s\n', ...
                pi, param, meanPairTune(pi), medianPairTune(pi), timing.tune.model(pi) + timing.tune.sample(pi));

            [~,minI] = min(meanPairTune);
            if(minI == pi)
                bestTuneModel = model;
            end

            if(meanPairTune(pi) >= 2*meanPairTune(1))
                fprintf('    Stopping tuning early since meanPairTune(cur) > 2*meanPairTune(1)\n\n');
                break;
            end
        end

        %% Display results
        fprintf('    %10s, %10s, %10s\n','paramVec', 'meanPair', 'medianPair');
        fprintf('    %10.5g, %10.5g, %10.5g\n', [paramVec;meanPairTune;medianPairTune]);
        [~,minI] = min(meanPairTune);
        fParam = paramVec(minI);
        timing.tune.total = toc(ts);
        fprintf('\n    Selected parameter = %d with mean pair value of = %g\n',fParam,meanPairTune(minI));
    end
    fprintf('  << Finished hyperparameter tuning >>\n\n');
    
    %% Use final param value to train 
    ts = tic;
    fModel = trainFunc( XtTrain, fParam, bestTuneModel );
    timing.train = toc(ts);
    
    %% Sample from final model
    ts = tic;
    fXtSample = sampleFunc( fModel, nSamples );
    timing.sample = toc(ts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple wrapper functions for sampling or model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XtSample = wrapperSampleSQR( PsiExt, nSamples, nGibbs, nInner, nBatches )
    XtSample = cell(nBatches,1);
    nSamplesVec = floor(nSamples/nBatches)*ones(1,nBatches);
    oneMore = nSamples - sum(nSamplesVec);
    nSamplesVec(1:oneMore) = nSamplesVec(1:oneMore) + 1;
    assert(sum(nSamplesVec) == nSamples, 'Batches not split properly');
    assert(all(nSamplesVec == floor(nSamples/nBatches) | nSamplesVec == floor(nSamples/nBatches)+1), 'Batches not split properly');

    parfor bi = 1:nBatches;
        nSamplesCur = nSamplesVec(bi);
        XtSample{bi} = mrfs.grm.univariate.Poisson.sampleSQR_Gibbs( PsiExt, nSamplesCur, nGibbs, nInner );
    end
    XtSample = cell2mat(XtSample);
end

function XtSample = wrapperSampleMixture( model, nSamples )
    % Extract parameters from model
    poissMean = model.poissMean;
    pVec = model.pVec;
    [k, p] = size(poissMean);
    
    % Sample from multinomial
    nFromCluster = mnrnd(nSamples, pVec);
    XtSample = NaN(nSamples, p);
    iCur = 0;
    
    % Sample from Poisson distributions
    for j = 1:k
        curN = nFromCluster(j);
        XtSample((iCur+1):(iCur+curN),:) = poissrnd(repmat(poissMean(j,:),curN,1));
        iCur = iCur + curN;
    end
end

function model = wrapperFlpgm( Xt, param, modLParam )
    % Train negative binomial on counts
    Lvec = full(sum(Xt, 2));
    dispersion = var(Lvec)/mean(Lvec);
    if(dispersion > 1)
        params = nbinfit(Lvec);
    else
        params = [NaN, NaN];
    end
    poissMean = mean(Lvec);

    % Train LPMRF model
    lpmrf = mrfs.models.LPMRF();
    beta = 1e-4;
    learner = mrfs.learners.PoissonRegModL( param, beta, modLParam );
    lpmrf.train( Xt, learner );
    
    warning('off','all');
    model = struct(lpmrf);
    warning('on','all');
    model.lpmrfObj = lpmrf;
    model.negbinR = params(1);
    model.negbinP = params(2);
    model.poissMean = poissMean;
    model.beta = beta;
    model.modLParam = modLParam;
    model.dispersion = dispersion;
end

function XtSample = wrapperSamplerFlgpm( model, nSamples, nAnneal, nGibbs, useNegBin )
    if(model.dispersion > 1 && useNegBin )
        LvecSample = nbinrnd(model.negbinR, model.negbinP, nSamples, 1);
    else
        LvecSample = poissrnd(model.poissMean, nSamples, 1);
    end
    
    p = length(model.thetaNode);
    curI = 0;
    XtSample = NaN(nSamples, p);
    for L = 0:max(LvecSample)
        nL = sum(LvecSample==L);
        sampler = mrfs.samplers.LPMRF_AISSampler(model.lpmrfObj, nAnneal, nGibbs );
        [tempSample, logW] = sampler.sampleL( nL, L, false);

        % Resample to get equally weighted samples
        w = exp(logW - max(logW));
        w = w./sum(w);
        resampleI = mnrnd(nL, w);
        for ii = 1:length(resampleI)
            idx = (curI+1):(curI+resampleI(ii));

            XtSample(idx,:) = repmat(tempSample(ii,:), resampleI(ii), 1);
            curI = curI+resampleI(ii);
        end
        %fprintf('Finished L = %d/%d\n', L, max(LvecSample));
    end
end

function model = wrapperPgm( Xt, lam )
    [n,p] = size(Xt);
    
    thetaNode = zeros(p,1);
    thetaEdge = zeros(p,p);
    parfor s = 1:p
        maxit = 1e6; thr = 1e-6;
        sel = false(p,1);
        sel(s) = true;
        [thetaNode(s), temp, obj, ind, iter] = EXP_GM_NeighborLearning('poisson', Xt(:,sel), Xt(:,~sel), lam, maxit, thr);
        thetaEdgeS = NaN(p,1);
        thetaEdgeS(~sel) = temp;
        thetaEdgeS(sel) = 0;
        thetaEdge(:,s) = thetaEdgeS;
        fprintf('s=%d, iter=%d, obj=%g, ind=%g\n', s, iter, obj(end), ind);
    end
    
    model.thetaNode = thetaNode;
    model.thetaEdge = (thetaEdge + thetaEdge')/2; % Symmetric and divide by 2 since
end
