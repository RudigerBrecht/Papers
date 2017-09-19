function udry=rbfExtrapolationLocalPolynomial(xwet,uwet,xdry,eps,nn,dx)

rscale=dx;
    
if isempty(xdry)
    udry=zeros(0,1);
    return;
end

kdTree=KDTreeSearcher(xwet);
idx=knnsearch(kdTree,xdry,'k',nn);

udry=zeros(nn,1);

for ii=1:length(xdry)
    
    o=ones(nn,1);
    xe=xwet(idx(ii,:));
    
    r=sqrt((o*xe'-xe*o').^2);
    
    phi=ones(nn+2);
    phi(nn+1:end,nn+1:end)=0;    
    phi(1:nn,1:nn)=sqrt(1+(eps*r).^2);
    phi(1:nn,nn+2)=xe';
    phi(nn+2,1:nn)=xe;

    w=phi\[uwet(idx(ii,:));0;0];

    rr=sqrt((xe-xdry(ii)).^2);
    phi=sqrt(1+(eps*rr/rscale).^2);

    udry(ii)=sum(w(1:nn).*phi)+w(nn+1)+w(nn+2)*xdry(ii);
   
end

end
