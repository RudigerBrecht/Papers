function [M]=computeAveragingMatrix(x,nn)

n=length(x);
M=zeros(n);

l=1;

av=sort(gaussianAverage(nn),'descend');

for ii=1:n
      
     Mx=zeros(1,n);
            
     xc=x(ii);

     % Find the nn nearest neighbors
     [~,idx]=sort(abs(xc-x));
     idx=idx(1:nn);
     % Use Gaussian averaging
     Mx(idx)=av; 
                                 
     M(l,:)=Mx(:)';           
     l=l+1;

end

M=sparse(M);

end

function gaussFilter=gaussianAverage(n)

sigma = 1;
size = n;
x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter);

end