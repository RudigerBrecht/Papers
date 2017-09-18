function lakeAtRest2d()

% integration time and timestepping
t=20;
dt=0.0015;


steps=round(t/dt);

n=50;
m=n;
x=linspace(0,1,n)';

y=x;


dx=x(2)-x(1); dy=y(2)-y(1);

[X,Y]=meshgrid(x,y);
X=X'; Y=Y';

ix=2:m-1; iy=2:n-1;

X(ix,iy)=X(ix,iy)+0.1*dx*randn(m-2,n-2); %+0.1*
Y(ix,iy)=Y(ix,iy)+0.1*dy*randn(m-2,n-2); %+0.1*




% Interpolation threshold
delta=0.05;%0.01;%0.00025;

% Gravitational constant
g=9.81;

% Bottom topography
b=zeros(n,n);
rm=0.4;
a=1.2;


for ii=1:n
    for jj=1:n
        r=sqrt((x(ii)-0.5)^2+(y(jj)-0.5)^2);        
%         r=sqrt((x(ii))^2+(y(jj))^2);
        if(r<0.4)
            b(ii,jj)=a*(exp(-0.5/(rm^2-r^2)))/(exp(-0.5/rm^2));
        end
    end
end


% Initial conditions
u=zeros(m,n); v=zeros(m,n);

h=ones(m,n);
  
h=h-b;
h0=h;


% Extrapolation epsilon
epsE=1;%5;%0.5;

% Number of nearest neighbors for local RBF extrapolation
nnE=40;%30

% Differentation matrix
eps=1; nn=25; % 0.1   9
[Dx,Dy,Lap]=rbfFD2dGA(X,Y,eps,nn,dx,dy);


% Compute averaging matrix
Ma=computeAveragingMatrix2d(X,Y,nn,dx,dy);


% L2 errors
normH=zeros(steps,1);
normH(1)=0;
normHu=normH;


% Bottom topography derivatives
bx=Dx*b(:); by=Dy*b(:);

% Viscosity
nu=0.003;

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
    hmean=Ma*h(:);    
   h2x=hmean.*(Dx*h(:));
   h2y=hmean.*(Dy*h(:));   
   
   
   % RK2 scheme
   
   % Euler forward step (half step)      
   uhn(wet)=u(wet).*h(wet)-...
       dt/2*(u2hx(wet)+uvhy(wet)+g*h2x(wet)+g*hmean(wet).*bx(wet)-nu*Del2u(wet));  
   vhn(wet)=v(wet).*h(wet)-...
       dt/2*(v2hy(wet)+uvhx(wet)+g*h2y(wet)+g*hmean(wet).*by(wet)-nu*Del2v(wet));  
   hn(wet)=h(wet)-dt/2*(uhx(wet)+vhy(wet));
   un(wet)=uhn(wet)./hn(wet); vn(wet)=vhn(wet)./hn(wet);
           
   
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
   hmean=Ma*hn(:);   
   h2x=hmean.*(Dx*hn(:));
   h2y=hmean.*(Dy*hn(:));
   
   % Euler forward step         
   uhn(wet)=u(wet).*h(wet)-...
       dt*(u2hx(wet)+uvhy(wet)+g*h2x(wet)+g*hmean(wet).*bx(wet)-nu*Del2un(wet));  
   vhn(wet)=v(wet).*h(wet)-...
       dt*(v2hy(wet)+uvhx(wet)+g*h2y(wet)+g*hmean(wet).*by(wet)-nu*Del2vn(wet));  
   hn(wet)=h(wet)-dt*(uhx(wet)+vhy(wet));
   un(wet)=uhn(wet)./hn(wet); vn(wet)=vhn(wet)./hn(wet);
        

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
   
   
   % Change in water height
    normH(ii)=norm(h(wet)-h0(wet),inf);
    normHu(ii)=norm(uhn(wet),inf);
    
   
end

plot(linspace(0,t,steps),normH, 'black');
xlabel('t')
ylabel('|| h-h_0 ||_\infty')

figure
plot(linspace(0,t,steps),normHu,'black');
xlabel('t')
ylabel('|| hu-(hu)_0 ||_\infty')

figure
plot(X(:),Y(:),'ko');
xlabel('x'); ylabel('y')
set(gca,'fontsize',14)
grid
title('Node distribution')



end