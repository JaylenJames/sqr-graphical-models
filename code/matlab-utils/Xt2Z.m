function Zemp = Xt2Z( Xt, maxX )
    %% Compute empirical probability matrix Z
    [n,p] = size(Xt);
    empHist = zeros(repmat(maxX+1,1,2));
    for i = 1:n
        ii = min(Xt(i,1)+1,size(empHist,1));
        jj = min(Xt(i,2)+1,size(empHist,2));
        empHist(ii,jj) = empHist(ii, jj) + 1;
    end
    Zemp = (empHist+eps)/sum(empHist(:));
end
