function [Ma]=computeAveragingMatrix2d(X,Y,nn,dx,dy)


[m,n]=size(X);
Ma=sparse(m*n,m*n);

% Inner domain
l=1;

kdTree = KDTreeSearcher([X(:),Y(:)]);

for jj=1:n
   for ii=1:m
       
       Mxy=zeros(m,n);

       xc=X(ii,jj); yc=Y(ii,jj);
       idx=knnsearch(kdTree,[xc,yc],'k',nn);
                      
       Xs=X(idx); Ys=Y(idx);  
       xs=Xs(:); ys=Ys(:);
           
       
       Mxy(idx)=exp(-(((xs-xc)/dx).^2+((ys-yc)/dy).^2));
       Mxy(idx)=Mxy(idx)./sum(Mxy(idx));
              
       Ma(l,:)=Mxy(:)';        
       l=l+1;
   end
end

end