function [Dx,Del2]=rbfFD1GApoly(x,eps,nn,dx)

if nargin<4
    dx=1;
end

n=length(x);
Dx=zeros(n);
Del2=zeros(n);

rscale=dx;

% Inner domain
l=1;
o=ones(nn,1);

for ii=1:n
      
     Lx=zeros(1,n); LDel2=zeros(1,n);
            
     xc=x(ii);

     % Find the nn nearest neighbors
     [~,idx]=sort(abs(xc-x));
     idx=idx(1:nn);
     
     xs=x(idx)';
         
     phi=zeros(nn+2);
     phi(nn+1,1:nn)=1; phi(1:nn,nn+1)=1;
     phi(nn+2,1:nn)=xs';  phi(1:nn,nn+2)=xs;
          
     r=sqrt((o*xs'-xs*o').^2);
     phi(1:nn,1:nn)=exp(-(eps*r/rscale).^2);
   
     rr=sqrt((xc*o-xs).^2);
     phix=-2*eps^2/rscale^2*(xc*o-xs).*exp(-(eps*rr/rscale).^2);
     phiDel2=eps^2/rscale^2*(4*(eps*rr/rscale).^2-4).*exp(-(eps*rr/rscale).^2);
           
     wx=phi\[phix;0;1];
     wDel2=phi\[phiDel2;0;0];
          
     Lx(idx)=wx(1:nn); LDel2(idx)=wDel2(1:nn);
                                 
     Dx(l,:)=Lx(:)';  Del2(l,:)=LDel2(:)';
     l=l+1;

end

Dx=sparse(Dx);

end