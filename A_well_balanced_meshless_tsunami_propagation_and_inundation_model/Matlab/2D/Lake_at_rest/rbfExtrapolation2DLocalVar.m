function varargout=rbfExtrapolation2DLocalVar(xwet,ywet,xdry,ydry,eps,nnE,dx,dy,varargin)

nVarargs=length(varargin);

if isempty(xdry)
    for kk=1:nVarargs
      varargout{kk}=zeros(0,1);
    end
end

rscale=sqrt(dx^2+dy^2);
    
kdTree=KDTreeSearcher([xwet,ywet]);
idx=knnsearch(kdTree,[xdry,ydry],'k',nnE);

% Preallocate output variables
for kk=1:nVarargs
  varargout{kk}=zeros(nnE,1);
end

rInfluence=3*sqrt(dx^2+dy^2);

for ii=1:length(xdry)
    
    o=ones(nnE,1);
    xe=xwet(idx(ii,:)); ye=ywet(idx(ii,:));
    
    phi=ones(nnE+1);
    phi(nnE+1,nnE+1)=0;

    r=sqrt((o*xe'-xe*o').^2+(o*ye'-ye*o').^2);
    phi(1:nnE,1:nnE)=sqrt(1+(eps*r/rscale).^2);
    rr=sqrt((xe-xdry(ii)).^2+(ye-ydry(ii)).^2);
    phie=sqrt(1+(eps*rr/rscale).^2);
    
    for kk=1:nVarargs
      
        if (min(rr)<rInfluence)
          w=phi\[varargin{kk}(idx(ii,:));0];
          varargout{kk}(ii)=sum(w(1:nnE).*phie)+w(nnE+1);
        else
          varargout{kk}(ii)=0;  
        end
    end
     
end

end
