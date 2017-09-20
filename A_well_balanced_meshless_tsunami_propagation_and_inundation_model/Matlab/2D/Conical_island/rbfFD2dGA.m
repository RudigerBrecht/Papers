function [Dx,Dy,Lap]=rbfFD2dGA(X,Y,eps,nn,dx,dy)

if nargin <5
    dx=1/sqrt(2); dy=1/sqrt(2);
end

rscale=sqrt(dx^2+dy^2);

[m,n]=size(X);
Dx=sparse(m*n,m*n); Dy=sparse(m*n,m*n); Lap=sparse(m*n,m*n);

% Inner domain
l=1;
o=ones(nn,1);

kdTree = KDTreeSearcher([X(:),Y(:)]);

for jj=1:n
   for ii=1:m
       
       Lx=sparse(m,n); Ly=sparse(m,n); LLap=sparse(m,n);
       
%        if (ii==1 || ii==m || jj==1 || jj==n)
                     
%        else    
       
           xc=X(ii,jj); yc=Y(ii,jj);
           idx=knnsearch(kdTree,[xc,yc],'k',nn);
           
           Xs=X(idx); Ys=Y(idx);
         
           xs=Xs(:); ys=Ys(:);
           
           r=sqrt((o*xs'-xs*o').^2+(o*ys'-ys*o').^2);
           phi=ones(nn+1); phi(nn+1,nn+1)=0;
           phi(1:nn,1:nn)=exp(-(eps*r/rscale).^2);

           rr=sqrt((xc*o-xs).^2+(yc*o-ys).^2);

           phix=-2*eps^2/rscale^2*(xc*o-xs).*exp(-(eps*rr/rscale).^2);
           phiy=-2*eps^2/rscale^2*(yc*o-ys).*exp(-(eps*rr/rscale).^2);
           
           % Laplacian
           phiLap=eps^2/rscale^2*(4*(eps*rr/rscale).^2-4).*exp(-(eps*rr/rscale).^2);
           
           % Laplacian^2
%            phiLap=eps^4/rscale^4*(16*(eps*rr/rscale).^4-...
%                64*(eps*rr/rscale).^2+32).*exp(-(eps*rr/rscale).^2);
           
           wx=phi\[phix;0]; wy=phi\[phiy;0]; wLap=phi\[phiLap;0];
           Lx(idx)=wx(1:nn); Ly(idx)=wy(1:nn); LLap(idx)=wLap(1:nn);
                                 
%        end
       
          Dx(l,:)=Lx(:)'; Dy(l,:)=Ly(:)';  Lap(l,:)=LLap(:)';
          l=l+1;
   end
end

% Dx=sparse(Dx); Dy=sparse(Dy); Lap=sparse(Lap);

end