PGM.path.neighborhood.mod <-
function(X,Y,nlams=10,startb=0,lambda=NULL)
{
	n = nrow(X); p = ncol(X);
	if(is.null(lambda)){
		lmax = lambdaMax(X)
		lambda = exp(seq(log(lmax),log(0.001*lmax),l=nlams))
	}
	fprintf <- function(...) cat(sprintf(...));
	#lams = exp(seq(log(lmax),log(.0001),l=nlams));
	#-if(nlams==1){lams = lmax};
	thr = 1e-8; maxit = 1e4;
	Xt = cbind(t(t(rep(1,n))),X);
	L = max(eigen(t(X)%*%X,only.values=TRUE)$values);
    L2 = max(t(Xt)%*%Y);
    
	alphas = 0; Bmat = matrix(0,p,nlams);
	if(sum(startb)==0){Bhat = matrix(runif(p+1),p+1,1); Bhat[1] = -log(mean(Y)); }else{Bhat=startb}
	for(i in 1:nlams){
		iter = 1; ind = 1;
		while(thr<ind & iter<maxit){
			oldb = Bhat;
            
			#tmp = Bhat - (t(Xt)%*%Y - t(Xt)%*%exp(-Xt%*%Bhat))/L;
            A = t(Xt)%*%Y;
            B = exp(-Xt%*%Bhat);
            C = t(Xt)%*%B;
            grad = (A-C)/L;
            
            tmp = Bhat - grad;
			#Bhat = matrix(sapply(tmp - lams[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
			Bhat = matrix(sapply(tmp - lambda[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
			ind = sum((Bhat - oldb)^2);
            fprintf('iter = %d, ind = %g, minBhat = %g, maxB = %g, maxC = %g, maxGrad = %g, L = %g, L2 = %g, lambda = %g\n',iter,ind,min(Bhat),max(B),max(C),max(grad),L,L2,lambda);
            show(Bhat);
			iter = iter + 1;
		}
		alphas[i] = -Bhat[1];
		Bmat[,i] = -Bhat[2:(p+1),drop=FALSE];
	}
  #return(list(alpha=alphas,Bmat=Bmat,lams=lams))
stop('After 1 variable')
  return(list(alpha=alphas,Bmat=Bmat,lambda=lambda))
}
