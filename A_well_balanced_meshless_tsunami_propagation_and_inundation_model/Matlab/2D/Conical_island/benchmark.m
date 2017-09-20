function F=benchmark

epsilon=0.04; % Case "A"
% epsilon=0.18; % Case "C" 


t=20;
dt=0.02;    


steps=round(t/dt);

% Grid
m=200; n=200;
x=linspace(0,25,m);
y=linspace(0,28,n);
dx=x(2)-x(1); dy=y(2)-y(1);

[X,Y]=meshgrid(x,y);
X=X'; Y=Y';

% position of gages;  
  g6x=find(abs(x-9.36)<=0.063); % 0.063 %0.043
  g6y=find(abs(y-13.8)<=0.1);
  
  g9x=find(abs(x-10.36)<=0.063);
  g9y=find(abs(y-13.8)<=0.1);

  g16x=find(abs(x-12.96)<=0.08); %0.08 % 0.05
  g16y=find(abs(y-11.22)<=0.1);

  g22x=find(abs(x-15.56)<=0.08); % 0.08 % 0.05
  g22y=find(abs(y-13.8)<=0.1);

% Interpolation threshold
delta=1/70;%1/50;

% Gravitational constant
g=9.81;

% Inner domain indices
ix=2:m-1; iy=2:n-1;

% Bottom topography
b=zeros(m,n);
for ii=1:m
    for jj=1:n
        xc=x(ii);
        yc=y(jj);
        r=sqrt((xc-12.96)^2+(yc-13.8)^2);
        if(r>3.6)
            
        elseif(r<1.1)
            b(x==xc,y==yc)=0.625;
        else
           b(x==xc,y==yc)=(r - 3.6)/(1.1 - 3.6)*0.625;
        end                 
    end
end

h0=0.32;
d=h0;
        
waveAmplitude=h0*epsilon;

% Initial conditions
u=zeros(m,n); v=zeros(m,n);

h=h0*ones(m,n);


H0=0.32; H=0.32*epsilon; T = 2.45;C=sqrt(g*(H0 + H));


h=h-b;

h0=h;

 gages=zeros(steps,4); % columns #6|#9|#16|#22
 gages(1,1)=h(g6x,g6y); 
 gages(1,2)=h(g9x,g9y);
 gages(1,3)=h(g16x,g16y); 
 gages(1,4)=h(g22x,g22y); 

% Extrapolation epsilon
epsE=0.7;%0.5;

% Number of nearest neighbors for local RBF extrapolation
nnE=15; %10

% Differentation matrix
eps=0.1; nn=9;
% [Dx,Dy,Lap]=rbfFD2dGA(X,Y,eps,nn,dx,dy);

% save('Dx.mat','Dx'); save('Dy.mat','Dy');


% % Compute averaging matrix
% Ma=computeAveragingMatrix2d(X,Y,nn,dx,dy);
% MDx=innerProductMW(Ma,Dx);
% MDy=innerProductMW(Ma,Dy);
% save('Ma.mat','Ma');  save('MDx.mat','MDx'); save('MDy.mat','MDy'); save('Lap.mat','Lap');
% 
% return

% Load precomputed matrices
load('Ma.mat'); load('MDx.mat'); load('MDy.mat'); load('Dx.mat'); load('Dy.mat'); load('Lap.mat');
% 

% Conserved quantities
M=zeros(steps,1);
M(1)=sum(sum(h*dx*dy));

normH=zeros(steps,1);
normH(1)=0;

% Bottom topography derivatives
bx=Dx*b(:); by=Dy*b(:);

% Viscosity
nu=0.01;

% Find initial dry points
dry0=h<=delta;

% Initialize array to record maximum flooding
hflooded=zeros(m,n);
% to save the maximal runup
hmax=h;
dry1=h>delta;
hmax(dry1)=0;
hmax=sparse(hmax);

% Time stepping
un=u; vn=v; hn=h; uhn=u.*h; vhn=v.*h;
loops=steps/10;
F(loops) = struct('cdata',[],'colormap',[]);
l=1;


tic
for ii=2:steps
   
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

   
   eta0=H*sech(sqrt(3*H/(4*(H0)^3))*C*(ii*dt - T)).^2;

   hn(1,:)= H0 + eta0;
   un(1,:)=C*eta0/(H0 + eta0);
   vn(1,:)=0;
   
   hn(m,:)=(dx*(h(m,:)+b(m,:))+dt*c(m,:).*(h(m-1,:)+b(m-1,:)))./(dx+dt*c(m,:))-b(m,:);
   un(m,:)=u(m,:)+g*(h(m-1,:)+b(m-1,:)-h(m,:)-b(m,:))/dx*dt;
   vn(m,:)=0;
   
   hn(:,1)=(dx*(h(:,1)+b(:,1))+dt*c(:,1).*(h(:,2)+b(:,2)))./(dx+dt*c(:,1))-b(:,1);
   un(:,1)=0;
   vn(:,1)=v(:,1)-g*(h(:,2)+b(:,2)-h(:,1)-b(:,1))/dx*dt;
   
   hn(:,n)=(dx*(h(:,n)+b(:,n))+dt*c(:,n).*(h(:,n-1)+b(:,n-1)))./(dx+dt*c(:,n))-b(:,n);
   un(:,n)=0;
   vn(:,n)=v(:,n)+g*(h(:,n-1)+b(:,n-1)-h(:,n)-b(:,n))/dx*dt;
   

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

   
   hn(m,:)=(dx*(h(m,:)+b(m,:))+dt*c(m,:).*(h(m-1,:)+b(m-1,:)))./(dx+dt*c(m,:))-b(m,:);
   un(m,:)=u(m,:)+g*(h(m-1,:)+b(m-1,:)-h(m,:)-b(m,:))/dx*dt;
   vn(m,:)=0;
   
   hn(:,1)=(dx*(h(:,1)+b(:,1))+dt*c(:,1).*(h(:,2)+b(:,2)))./(dx+dt*c(:,1))-b(:,1);
   un(:,1)=0;
   vn(:,1)=v(:,1)-g*(h(:,2)+b(:,2)-h(:,1)-b(:,1))/dx*dt;
   
   hn(:,n)=(dx*(h(:,n)+b(:,n))+dt*c(:,n).*(h(:,n-1)+b(:,n-1)))./(dx+dt*c(:,n))-b(:,n);
   un(:,n)=0;
   vn(:,n)=v(:,n)+g*(h(:,n-1)+b(:,n-1)-h(:,n)-b(:,n))/dx*dt;

   % Local RBF extrapolation
   [hn(dry),un(dry),vn(dry)]=rbfExtrapolation2DLocalVar(X(wet),Y(wet),...
       X(dry),Y(dry),epsE,nnE,dx,dy,hn(wet)+b(wet),un(wet),vn(wet));     
    hn(dry)=hn(dry)-b(dry);
   
   h=hn; u=un; v=vn;
   
   % Record max runup
   flooded=find(wet==dry0);
   for kk=1:length(flooded)
     hflooded(flooded(kk))=max(hflooded(flooded(kk)),hn(flooded(kk)));
   end
   
   %find max runup
   htemp=h+b;
   htemp(dry1)=0;
   htemp=sparse(htemp);
   [idx,idy] = find(hmax);   
   hmax(idx,idy)=max(hmax(idx,idy),htemp(idx,idy));
   
   % collect data for gages
   gages(ii,1)=h(g6x,g6y)+b(g6x,g6y);
   gages(ii,2)=h(g9x,g9y)+b(g9x,g9y);
   gages(ii,3)=h(g16x,g16y)+b(g16x,g16y);
   gages(ii,4)=h(g22x,g22y)+b(g22x,g22y);
   

if epsilon==0.18   
    if(ii*dt==10 || ii*dt==14) 
       figure
        surf(X,Y,b);
        colormap gray
        freezeColors;
        hold on
        colormap cool
        mesh(X,Y,h+b);
        view(-30,70)
    end
end
if epsilon==0.04
     if(ii*dt==12 || ii*dt==16) 
       figure
        surf(X,Y,b);
        colormap gray
        freezeColors;
        hold on
        colormap cool
        mesh(X,Y,h+b);
        view(-30,70)
    end 
end
    
end


 save('gages.mat','gages');
 save('testRunUp.mat','hmax');
 
 plotData(epsilon)
 

end

