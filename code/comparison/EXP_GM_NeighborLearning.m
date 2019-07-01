%type: 'ising' or 'poisson' or 'gaussian' or 'exponential'
%Ys: n x 1 target variable (We will estimate the set of neighborhoods of random variable Ys)
%Ysc: n x (p-1) remaining variables in V\setmius s
%lam: regularization parameter
%maxit: maximum number of iteration of proximal graident algorithm
%thr: algorithm will stop if the 'change' of estimator after one iteration is less than thr 

function[Ahat,Bhat,obj,ind,iter] = EXP_GM_NeighborLearning(type,Ys,Ysc,lam,maxit,thr)
    %% ALGORITHM PARAMETERS
%     thr = 1e-6;
%     maxit = 1e6;
    beta = 0.5;
    minA = -0.001;
    
    %% PROBLEM SETUP
    [n,p] = size(Ysc);
    [n_temp,temp2] = size(Ys);

    assert(n==n_temp);
    assert(temp2==1);

    if(strcmp(type, 'ising')||strcmp(type, 'poisson')||strcmp(type, 'gaussian')||strcmp(type, 'exponential'))
    else
        fprintf(1, 'ERROR:distribution type should be one of ''isng'', ''poisson'', ''gaussian'' or ''exponential'' \n'); 
        return
    end
    
    yy = (repmat(Ys, [1 p]).*Ysc);
    
    %Bhat = rand(p,1);
    Bhat = zeros(p,1);
    if(strcmp(type,'poisson')||strcmp(type,'exponential'))
        Bhat = -1*Bhat;
    end
    %Ahat = rand();
    Ahat = log(max(mean(Ys),1e-5)); % Just to ensure non -inf value
    
    if(strcmp(type, 'exponential')) 
        Ahat = -Ahat;
    end
    
    iter = 1;
    ind = 1;
    obj = inf;

    while thr<ind & iter<maxit
        Bold = Bhat;
        Aold = Ahat;
        
        [D_Agrad D_Bgrad] = computeDgrad(type,Ysc,Ahat,Bhat);
        Bgrad = (-(1/n)*sum(yy,1) + (1/n)*D_Bgrad)';
        Agrad = -(1/n)*sum(Ys) + (1/n)*D_Agrad;
        
        t = 1;
  
        Btmp = Bold - t*Bgrad;
        Bhat = sign(Btmp).*max(abs(Btmp)-lam*t,0);
        
        if(strcmp(type,'poisson')||strcmp(type,'exponential')) %should be negative
            Bhat = min(Bhat,0);
        end

        G_b = (Bold-Bhat)/t;

        Ahat = Aold - t*Agrad;

        if(strcmp(type, 'exponential')) 
            Ahat = min(Ahat,minA);
        end
        
        G_a = (Aold-Ahat)/t;

        g_x =  -(1/n)*sum(yy*Bold) -(1/n)*sum(Ys)*Aold + (1/n)*computeD(type,Ysc,Aold,Bold);    
        
        while (-(1/n)*sum(yy*Bhat) -(1/n)*sum(Ys)*Ahat' + (1/n)*computeD(type,Ysc,Ahat,Bhat) ...
                    > g_x - t*(Bgrad(:)'*G_b(:)+Agrad*G_a') + t/2*(G_b(:)'*G_b(:)+G_a*G_a') ) 
                
            t = beta*t;

            Btmp = Bold - t*Bgrad;    
            Bhat = sign(Btmp).*max(abs(Btmp)-lam*t,0);
            
            if(strcmp(type,'poisson')||strcmp(type,'exponential'))
                Bhat = min(Bhat,0);
            end
        
            G_b = (Bold-Bhat)/t;

            Ahat = Aold - t*Agrad;
            
            if(strcmp(type, 'exponential'))
                Ahat = min(Ahat,minA);
            end
            
            G_a = (Aold-Ahat)/t;
        end

        iter = iter + 1;
        cur_obj = -(1/n)*sum(yy*Bhat) -(1/n)*sum(Ys)*Ahat + (1/n)*computeD(type,Ysc,Ahat,Bhat) + lam*sum(abs(Bhat));
        
        obj = [obj; cur_obj];
        ind = norm([Ahat(:);Bhat(:)] - [Aold(:);Bold(:)],2)/norm([Aold(:);Bold(:)],2);
    end
%     likelihood = -(1/n)*sum(yy*Bhat) -(1/n)*sum(Ys)*Ahat + (1/n)*computeD(type,Ysc,Ahat,Bhat);
end

function [Agrad Bgrad] = computeDgrad(type,Ysc,Ahat,Bhat)
    [n,p] = size(Ysc);
%     [n_temp,q] = size(X);

    eta = Ahat+Ysc*Bhat;
    
    if(strcmp(type, 'ising'))
        temp = exp( 2*eta );
        temp_Dgrad = (temp-1)./(temp+1) ;
    elseif(strcmp(type, 'poisson'))
        temp_Dgrad = exp( eta );
    elseif(strcmp(type, 'gaussian'))
        temp_Dgrad = eta;
    elseif(strcmp(type, 'exponential'))
        temp_Dgrad = -1./(eta);
    end
    
    Bgrad = sum(repmat(temp_Dgrad,[1 p]).*Ysc,1);
    Agrad = sum( temp_Dgrad);

end

function dVal = computeD(type,Ysc,Ahat,Bhat)
    eta = Ahat + Ysc*Bhat;
    if(strcmp(type, 'ising'))
        dVal = sum( log(exp(-eta)+exp(eta)) );
    elseif(strcmp(type, 'poisson'))
        dVal = sum( exp(eta) );
    elseif(strcmp(type, 'gaussian'))
        dVal = 0.5*sum(eta.^2);
    elseif(strcmp(type, 'exponential'))
        dVal = -sum(log(-eta));
    end
end
