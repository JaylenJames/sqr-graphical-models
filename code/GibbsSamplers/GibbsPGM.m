% Gibbs sampler for Poisson graphical model
function [X] = GibbsPGM(n,p,alpha,Th,maxit)
Theta=Th;

X=poissrnd(1,[n,p]);
dims = 1:size(X,2);
  iter = 1;
  while iter < maxit 
    for s = 1:p
        sIx = dims~=s;
        mu = exp(alpha(s) + X(:,sIx)*Theta(sIx,s));
        X(:,s) = poissrnd(mu);
    end
      iter = iter + 1;
  end
 end