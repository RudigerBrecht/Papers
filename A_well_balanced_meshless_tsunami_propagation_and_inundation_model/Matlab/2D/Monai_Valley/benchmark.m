function F=benchmark

tFinal=25;
dt=0.01;    

steps=round(tFinal/dt);

% Load bathymetry
% load topoMonaiValley.mat
load topoMonaiValleyHighRes.mat
fac=1;
% h = 1/9*ones(3);
% b=filter2(h,b);
X=X(1:fac:end,1:fac:end);
Y=Y(1:fac:end,1:fac:end);
b=b(1:fac:end,1:fac:end);
X=X'; Y=Y'; b=b';
dx=X(2,1)-X(1,1);
dy=Y(1,2)-Y(1,1);
[m,n]=size(b);

% Load incident wave
load incidentWave.mat

% Load gage data
load gages.mat;

% plot(tInput,hInput,'x')
% return

% surf(X,Y,b);
% return

% Interpolation threshold
delta=1/125;
% delta=1/1000;

% Gravitational constant
g=9.81;

% Inner domain indices
ix=2:m-1; iy=2:n-1;

% Initial conditions
u=zeros(m,n); v=zeros(m,n); h=zeros(m,n);
h=h-b;

% Load initial conditions
% load ICs.mat;
% load step1400.mat
% load Height19.5.mat

% Find gage positions
indGage1=find(sqrt((X(:)-4.521).^2+(Y(:)-1.196).^2)==min(sqrt((X(:)-4.521).^2+(Y(:)-1.196).^2)));
indGage2=find(sqrt((X(:)-4.521).^2+(Y(:)-1.696).^2)==min(sqrt((X(:)-4.521).^2+(Y(:)-1.696).^2)));
indGage3=find(sqrt((X(:)-4.521).^2+(Y(:)-2.196).^2)==min(sqrt((X(:)-4.521).^2+(Y(:)-2.196).^2)));

hGage1=zeros(steps,1); hGage2=zeros(steps,1); hGage3=zeros(steps,1);
hGage1(1)=h(indGage1)+b(indGage1); hGage2(1)=h(indGage2)+b(indGage2); hGage3(1)=h(indGage3)+b(indGage3);

% h=hn;
% surf(X,Y,b);
% colormap gray
% freezeColors;
% hold on
% colormap cool
% mesh(X,Y,h+b);
% drawnow;
% hold off;
% return

% Extrapolation epsilon
epsE=5;

% Number of nearest neighbors for local RBF extrapolation
nnE=20;

% Differentation matrix
eps=0.1; nn=9;
% [Dx,Dy,Lap]=rbfFD2dGA(X,Y,eps,nn,dx,dy);

% Compute averaging matrix
% Ma=computeAveragingMatrix2d(X,Y,nn,dx,dy);
% MDx=innerProductMW(Ma,Dx); MDy=innerProductMW(Ma,Dy);
% save('MDx.mat','MDx'); save('MDy.mat','MDy');
% save('Ma.mat','Ma'); save('Lap.mat','Lap'); save('Dx.mat','Dx'); save('Dy.mat','Dy');

% load('MDx.mat'); load('MDy.mat'); 
load('Ma.mat'); load('Dx.mat'); load('Dy.mat'); load('Lap.mat');

% return

% h(1,:)=h(1,:)+0.15;

% Conserved quantities
M=zeros(steps,1);
wet=h>delta;
M(1)=sum(sum(h(wet)*dx*dy));

% Bottom topography derivatives
bx=Dx*b(:); by=Dy*b(:);

% Viscosity
nu=0.0025; % works
% nu=0.004; % works as well
% nu=0.001;

% Time stepping
un=u; vn=v; hn=h; uhn=u.*h; vhn=v.*h;
loops=steps/10;
F(loops) = struct('cdata',[],'colormap',[]);
l=1; k=2;
for ii=2:steps
   
   tcurr=(ii-1)*dt;
   if (k<=length(tInput))
   if (abs(tcurr-tInput(k))<1e-10)
       h(1,:)=hInput(k)-b(1,:);
       k=k+1;
   end
   end
   
   fprintf('Compute %ith step out of %i steps.\n',ii,steps);
   
   wet=h>delta; dry=h<=delta; 
            
   % Derivatives for conservative form
   u2hx=Dx*(u(:).^2.*h(:));  uhx=Dx*(h(:).*u(:)); 
   vhy=Dy*(h(:).*v(:));  uvhy=Dy*(h(:).*u(:).*v(:)); 
   uvhx=Dx*(h(:).*u(:).*v(:));  v2hy=Dy*(v(:).^2.*h(:));
   Del2u=Lap*(h(:).*u(:));  Del2v=Lap*(h(:).*v(:));
   
   % Averaging matrix
      
   % Mean value for h
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
              
   % Phase speed for outflow boundary conditions
   c=0*h; c(wet)=sqrt(g*h(wet));

   % Boundary conditions (Outflow)
   hn(1,:)=(dx*(h(1,:)+b(1,:))+dt*c(1,:).*(h(2,:)+b(2,:)))./(dx+dt*c(1,:))-b(1,:);
   un(1,:)=u(1,:)-g*(h(2,:)+b(2,:)-h(1,:)-b(1,:))/dx*dt;
   vn(1,:)=0;
   
   hn(m,:)=hn(m-1,:)+b(m-1,:)-b(m,:);
   un(m,:)=-un(m-1,:); vn(m,:)=vn(m-1,:);
   
   hn(:,1)=(hn(:,2)+b(:,2))-b(:,1);
   vn(:,1)=-vn(:,2); un(:,1)=un(:,2);
   
   hn(:,n)=(hn(:,n-1)+b(:,n-1))-b(:,n);
   vn(:,n)=-vn(:,n-1); un(:,n)=un(:,n-1);
      
   % Local RBF extrapolation
   [hn(dry),un(dry),vn(dry)]=rbfExtrapolation2DLocalVar(X(wet),Y(wet),...
       X(dry),Y(dry),epsE,nnE,dx,dy,hn(wet)+b(wet),un(wet),vn(wet));     
    hn(dry)=hn(dry)-b(dry);
   
   % Runge Kutta step (2nd order)
   
   % Derivatives for conservative form
   u2hx=Dx*(un(:).^2.*hn(:));  uhx=Dx*(hn(:).*un(:)); 
   vhy=Dy*(hn(:).*vn(:));  uvhy=Dy*(hn(:).*un(:).*vn(:)); 
   uvhx=Dx*(hn(:).*un(:).*vn(:));  v2hy=Dy*(vn(:).^2.*hn(:));
   Del2un=Lap*(hn(:).*un(:));  Del2vn=Lap*(hn(:).*vn(:));
   
  % Mean value for h
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
           
   hn(1,:)=(dx*(h(1,:)+b(1,:))+dt*c(1,:).*(h(2,:)+b(2,:)))./(dx+dt*c(1,:))-b(1,:);
   un(1,:)=u(1,:)-g*(h(2,:)+b(2,:)-h(1,:)-b(1,:))/dx*dt;
   vn(1,:)=0;
   
   hn(m,:)=hn(m-1,:)+b(m-1,:)-b(m,:);
   un(m,:)=-un(m-1,:); vn(m,:)=vn(m-1,:);
   
   hn(:,1)=(hn(:,2)+b(:,2))-b(:,1);
   vn(:,1)=-vn(:,2); un(:,1)=un(:,2);
   
   hn(:,n)=(hn(:,n-1)+b(:,n-1))-b(:,n);
   vn(:,n)=-vn(:,n-1); un(:,n)=un(:,n-1);

   % Local RBF extrapolation
   [hn(dry),un(dry),vn(dry)]=rbfExtrapolation2DLocalVar(X(wet),Y(wet),...
       X(dry),Y(dry),epsE,nnE,dx,dy,hn(wet)+b(wet),un(wet),vn(wet));     
    hn(dry)=hn(dry)-b(dry);
   
   h=hn; u=un; v=vn;
      
      
   % Mass
%    M(ii)=sum(sum(h(wet)*dx*dy));
      
   % Gage heights
   hGage1(ii)=h(indGage1)+b(indGage1); 
   hGage2(ii)=h(indGage2)+b(indGage2); 
   hGage3(ii)=h(indGage3)+b(indGage3);
   
end


t=linspace(0,tFinal,length(hGage1));

figure
subplot(311); plot(t,hGage1,'k');
hold on
subplot(311); plot(tGages(tGages<=tFinal),gage1Ref(tGages<=tFinal)/10,'k--');
ylabel('Gauge1'); grid
set(gca,'fontsize',14);
subplot(312); plot(t,hGage2,'k');
hold on
subplot(312); plot(tGages(tGages<=tFinal),gage2Ref(tGages<=tFinal)/10,'k--');
ylabel('Gauge2'); grid
set(gca,'fontsize',14);
subplot(313); plot(t,hGage3,'k');
set(gca,'fontsize',14);
hold on
subplot(313); plot(tGages(tGages<=tFinal),gage3Ref(tGages<=tFinal)/10,'k--');
ylabel('Gauge3'); grid
xlabel('t')

set(gca,'fontsize',14);

end

