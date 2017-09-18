function lakeAtRest2dtest(ep)

if nargin <1

end

t=20;%20;%20
dt=0.0015;%0.0015;%0.0025;


tic
steps=round(t/dt);

n=50;
m=n;
x=linspace(0,2,n)';
y=linspace(0,1,m)';
%  x=linspace(0,5,n)';
% y=x;


dx=x(2)-x(1); dy=y(2)-y(1);

[X,Y]=meshgrid(x,y);
X=X'; Y=Y';

ix=2:m-1; iy=2:n-1;

X(ix,iy)=X(ix,iy)+0.1*dx*randn(m-2,n-2); %+0.1*
Y(ix,iy)=Y(ix,iy)+0.1*dy*randn(m-2,n-2); %+0.1*

% save('X.mat','X');
% save('Y.mat','Y');

% load('X.mat');
% load('Y.mat');

% plot(X(:),Y(:),'ko');
% xlabel('x'); ylabel('y')
% set(gca,'fontsize',14)
% grid
% return

% Interpolation threshold
delta=0.01;%0.01;%0.00025;

% Gravitational constant
g=9.81;

% Bottom topography
bs=zeros(n,n);
rm=0.4;
a=1.2;


for ii=1:n
    for jj=1:n
        bs(ii,jj)=0.8*exp(-5*(x(ii)-0.9)^2-50*(y(jj)-0.5)^2);
    end
end

bn=zeros(n,1);
a=[0.1,0.2,0.3];
p=[1.6,3.2,0.5]';
for ii=1:n
    for jj=1:n
        bn(ii,jj)=bs(ii,jj)+0.5*(sum(a*sin(16*[1;2;3]*pi*X(ii,jj)+p)));%+sum(a*sin(16*[1;2;3]*pi*Y(ii,jj)+p)));
    end
end



%    b=zeros(m,n);

% Initial conditions
u=zeros(m,n); v=zeros(m,n);


h=2*ones(m,n);

 
 b=bs;  
% b=bn;

b=4*b;

h=h-b;

h0=h;

 
% hn=h;
% surf(X,Y,b);
% colormap gray
% freezeColors;
% hold on
% colormap winter
% temp=hn+b; temp(hn<delta)=nan;
% hp=surf(X,Y,temp);
% % axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(b(:)) max(b(:))])
% xlabel('x'); ylabel('y'), zlabel('h+b'); set(gca,'fontsize',14);
% % view(218,40);
% shading interp;
% lightangle(90,0);
% hp.FaceLighting = 'gouraud';
% hp.AmbientStrength = 0.3;
% hp.DiffuseStrength = 0.8;
% hp.SpecularStrength = 0.9;
% hp.SpecularExponent = 25;
% hp.BackFaceLighting = 'unlit';
% return

% Extrapolation epsilon
epsE=0.5;%5;%0.5;

% Number of nearest neighbors for local RBF extrapolation
nnE=10;%20;

% Differentation matrix
eps=10; nn=25; % 0.1   9
[Dx,Dy,Lap]=rbfFD2dGA(X,Y,eps,nn,dx,dy);

% Set animation configurations
%  f=figure;
%  set(f,'Visible', 'off');%'animating in background
%  a=gca;
%  %set(a,'DataAspectRatio',[1 1 1])
%  hold on
%  writerObj = VideoWriter('tsunami');
%  writerObj.Quality = 100;
%  writerObj.FrameRate = 30;%frames per sec
%  open(writerObj); 



% % Compute averaging matrix
Ma=computeAveragingMatrix2d(X,Y,nn,dx,dy);
MDx=innerProductMW(Ma,Dx);
MDy=innerProductMW(Ma,Dy);
% save('Ma.mat','Ma');  save('MDx.mat','MDx'); save('MDy.mat','MDy'); save('Lap.mat','Lap');
% save('Dx.mat','Dx'); save('Dy.mat','Dy');
% 
% return

% Load precomputed matrices
% load('Ma.mat'); load('MDx.mat'); load('MDy.mat'); load('Dx.mat'); load('Dy.mat'); load('Lap.mat');
% 

% Conserved quantities
% M=zeros(steps,1);
% M(1)=sum(sum(h*dx*dy));
% 
normH=zeros(steps,1);
normH(1)=0;
normHu=normH;

% Bottom topography derivatives
bx=Dx*b(:); by=Dy*b(:);

% Viscosity
nu=0.002;%0.001;

% Time stepping
un=u; vn=v; hn=h; uhn=u.*h; vhn=v.*h;
for ii=2:steps
   
    fprintf('Compute %ith step out of %i steps.\n',ii,steps);
   
    wet=h>delta; dry=h<=delta; 
        
   % Derivatives for conservative form
   u2hx=Dx*(u(:).^2.*h(:));  uhx=Dx*(h(:).*u(:)); 
   vhy=Dy*(h(:).*v(:));  uvhy=Dy*(h(:).*u(:).*v(:)); 
   uvhx=Dx*(h(:).*u(:).*v(:));  v2hy=Dy*(v(:).^2.*h(:));
   Del2u=Lap*((h(:)+b(:)).*u(:));  Del2v=Lap*((h(:)+b(:)).*v(:));
   
   % Averaging matrix
   h2x=computeDh2(MDx,h);  h2y=computeDh2(MDy,h);      
   h2x=reshape(h2x,m,n); h2y=reshape(h2y,m,n);
   
   hmean=reshape(Ma*h(:),m,n);
      
   % RK2 scheme
   
   % Euler forward step (half step)      
   uhn(wet)=u(wet).*h(wet)-...
       dt/2*(u2hx(wet)+uvhy(wet)+g*h2x(wet)+g*hmean(wet).*bx(wet)-nu*Del2u(wet));  
   vhn(wet)=v(wet).*h(wet)-...
       dt/2*(v2hy(wet)+uvhx(wet)+g*h2y(wet)+g*hmean(wet).*by(wet)-nu*Del2v(wet));  
   hn(wet)=h(wet)-dt/2*(uhx(wet)+vhy(wet));
   un(wet)=uhn(wet)./hn(wet); vn(wet)=vhn(wet)./hn(wet);
           
   % Global RBF extrapolation
%    hn(dry)=rbfExtrapolation2D(X(wet),Y(wet),hn(wet)+b(wet),X(dry),Y(dry),epsE)-b(dry);
%    un(dry)=rbfExtrapolation2D(X(wet),Y(wet),un(wet),X(dry),Y(dry),epsE);
%    vn(dry)=rbfExtrapolation2D(X(wet),Y(wet),vn(wet),X(dry),Y(dry),epsE); 
   
     % Boundary conditions (Reflective)
    un(1,:)=-un(2,:); un(n,:)=-un(n-1,:); vn(:,1)=-vn(:,2); vn(:,n)=-vn(:,n-1); 
    un(:,1)=un(:,2); un(:,n)=un(:,n-1); vn(1,:)=vn(2,:); vn(n,:)=vn(n-1,n);
    hn(1,:)=hn(2,:)+b(2,:)-b(1,:);  hn(n,:)=hn(n-1,:)+b(n-1,:)-b(n,:);
    hn(:,1)=hn(:,2)+b(:,2)-b(:,1);  hn(:,n)=hn(:,n-1)+b(:,n-1)-b(:,n);
    
   % Local RBF extrapolation
      [hn(dry),un(dry),vn(dry)]=rbfExtrapolation2DLocalVar(X(wet),Y(wet),...
          X(dry),Y(dry),epsE,nnE,dx,dy,hn(wet)+b(wet),un(wet),vn(wet));     
       hn(dry)=hn(dry)-b(dry);
   
   % Runge Kutta step (2nd order)
   
   % Derivatives for conservative form
   u2hx=Dx*(un(:).^2.*hn(:));  uhx=Dx*(hn(:).*un(:)); 
   vhy=Dy*(hn(:).*vn(:));  uvhy=Dy*(hn(:).*un(:).*vn(:)); 
   uvhx=Dx*(hn(:).*un(:).*vn(:));  v2hy=Dy*(vn(:).^2.*hn(:));
   Del2un=Lap*((hn(:)+b(:)).*un(:));  Del2vn=Lap*((hn(:)+b(:)).*vn(:));
   
   % Averaging matrix
   h2x=computeDh2(MDx,hn);  h2y=computeDh2(MDy,hn);      
   h2x=reshape(h2x,m,n); h2y=reshape(h2y,m,n);
   
   hmean=reshape(Ma*hn(:),m,n);
      
   % Euler forward step         
   uhn(wet)=u(wet).*h(wet)-...
       dt*(u2hx(wet)+uvhy(wet)+g*h2x(wet)+g*hmean(wet).*bx(wet)-nu*Del2un(wet));  
   vhn(wet)=v(wet).*h(wet)-...
       dt*(v2hy(wet)+uvhx(wet)+g*h2y(wet)+g*hmean(wet).*by(wet)-nu*Del2vn(wet));  
   hn(wet)=h(wet)-dt*(uhx(wet)+vhy(wet));
   un(wet)=uhn(wet)./hn(wet); vn(wet)=vhn(wet)./hn(wet);
        
   % Global RBF extrapolation
%     hn(dry)=rbfExtrapolation2D(X(wet),Y(wet),hn(wet)+b(wet),X(dry),Y(dry),epsE)-b(dry);
%     un(dry)=rbfExtrapolation2D(X(wet),Y(wet),un(wet),X(dry),Y(dry),epsE);
%     vn(dry)=rbfExtrapolation2D(X(wet),Y(wet),vn(wet),X(dry),Y(dry),epsE); 


     % Boundary conditions (Reflective)
    un(1,:)=-un(2,:); un(n,:)=-un(n-1,:); vn(:,1)=-vn(:,2); vn(:,n)=-vn(:,n-1); 
    un(:,1)=un(:,2); un(:,n)=un(:,n-1); vn(1,:)=vn(2,:); vn(n,:)=vn(n-1,n);
    hn(1,:)=hn(2,:)+b(2,:)-b(1,:);  hn(n,:)=hn(n-1,:)+b(n-1,:)-b(n,:);
    hn(:,1)=hn(:,2)+b(:,2)-b(:,1);  hn(:,n)=hn(:,n-1)+b(:,n-1)-b(:,n);

 
  
   % Local RBF extrapolation
      [hn(dry),un(dry),vn(dry)]=rbfExtrapolation2DLocalVar(X(wet),Y(wet),...
          X(dry),Y(dry),epsE,nnE,dx,dy,hn(wet)+b(wet),un(wet),vn(wet));     
       hn(dry)=hn(dry)-b(dry);
   
   h=hn; u=un; v=vn;
   

   
%         cla(a);   
% % %       
%         surf(X,Y,b);
%         colormap copper
%         freezeColors;
%         hold on
%         colormap jet
%         surf(X,Y,h+b);
% % %          view(-58,69)
% % %          view(0,30)
% %         view(-48,59)
%         drawnow;
%         hold off;
%          
%          frame = getframe(f);
%          writeVideo(writerObj,frame); 
       
%    M(ii)=sum(sum(h*dx*dy));
   
   % Change in water height
    normH(ii)=norm(h(wet)-h0(wet),inf);
    normHu(ii)=norm(uhn(wet),inf);
    
    plot(normH(2:ii))
    drawnow
   
end

toc
save('normH.mat','normH')
save('normHu.mat','normHu')

plot(linspace(0,t,steps),normH, 'black');
xlabel('t')
ylabel('|| h-h_0 ||_\infty')
figure
plot(linspace(0,t,steps),normHu,'black');
xlabel('t')
ylabel('|| hu-(hu)_0 ||_\infty')

return


end

function h2x=computeDh2(MW,h)

n=length(MW);

h2x=zeros(n,1);

for ii=1:n
%     h2x(ii)=h(:)'*((M(ii,:)'*W(ii,:)))*h(:);
    h2x(ii)=h(:)'*(MW{ii}*h(:));
end

end

function MW=innerProductMW(M,W)

n=length(M);

MW=cell(n,1);

for ii=1:n
   MW{ii}=M(ii,:)'*W(ii,:);
end


end