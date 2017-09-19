function benchmark2(n)

rng(1);
close all

t=3000;

if nargin <1
    n=200;
    dt=1;
end
   
if(n==200)
    dt=1;
end    
if(n==400)
    dt=1/2;
end
if(n==800)
    dt=1/4;
end
if(n==1600)
    dt=1/8;
end

    
steps=round(t/dt);

g=9.81;
L=5000;
x=linspace(-L,L,n)';
dx=x(2)-x(1);

% Threshold for dry areas
delta=1/100;

% Bottom topography
h0=10; a=3000;
b=h0*(x/a).^2;

% Initial conditions
omega=sqrt(2*g*h0)/a; B=5; t0=0;
u=B*a*omega/sqrt(2*h0*g)*sin(omega*t0)*ones(size(x));
h=h0-B^2/4/g*(1+cos(2*omega*t0))-B*x/2/a*sqrt(8*h0/g)*cos(omega*t0);
h=h-b;


% Differentation matrix
eps=0.1; nn=3;
[Dx,Del2]=rbfFD1GApoly(x',eps,nn,dx);
 

% Conserved quantities
M=zeros(steps,1);
wet=h>delta;
M(1)=sum(h(wet)*dx);

Herror=zeros(steps-1,1);
Huerror=zeros(steps-1,1);

% Time stepping
un=u; hn=h; uhn=u.*h;

nu=-0.005; epsE=1;

Ma=computeAveragingMatrix(x,nn);
% MDx=innerProductMW(Ma,Dx);

for ii=2:steps
    
   fprintf('Computing %ith step out of %i steps.\n',ii,steps);
   
   hana=h0-B^2/(4*g)*(1+cos(2*omega*(ii-1)*dt))-(B/(2*a))*x*sqrt(8*h0/g)*cos(omega*(ii-1)*dt)-b;
   uana=B*a*omega/(sqrt(2*h0*g))*sin(omega*(ii-1)*dt)*ones(size(x));
   
   wet=h>delta; dry=h<=delta;  
   
   % Derivatives for conservative form
   u2hx=Dx*(u.^2.*h); % h2x=0.5*Dx*(h.^2);
   uhx=Dx*(h.*u); bx=Dx*b; %hx=Dx*h;
   Del2u=Del2*(h.*u);
      
   % Averaging matrix
%    h2x=computeDxh2(MDx,h); 
           
   % Euler forward step
   hmean=Ma*h;
   h2x=hmean.*(Dx*h);
      
   uhn(wet)=u(wet).*h(wet)-dt*(u2hx(wet)+g*h2x(wet)+g*hmean(wet).*bx(wet)+nu*Del2u(wet));   
   hn(wet)=h(wet)-dt*uhx(wet);
   un(wet)=uhn(wet)./hn(wet);
         
   hn=extrapolatingBC(hn+b,x,dry,wet,epsE,nn,dx)-b;
   un=extrapolatingBC(un,x,dry,wet,epsE,nn,dx);
        
   % Boundary conditions (Reflective)
   un(1)=-un(2); un(n)=-un(n-1);
   hn(1)=hn(2)+b(2)-b(1);  hn(n)=hn(n-1)+b(n-1)-b(n);
      
   % Heun's method
                            
   % Derivatives for conservative form
   u2hxn=Dx*(un.^2.*hn); %h2xn=0.5*Dx*(hn.^2); 
   uhxn=Dx*(hn.*un); %hxn=Dx*hn;
   Del2un=Del2*(hn.*un);
              
   % Averaging matrix
%    h2xn=computeDxh2(MDx,hn);
   hmeann=Ma*hn;   
   h2xn=hmeann.*(Dx*hn);
       
   uhn(wet)=u(wet).*h(wet)-dt/2*(u2hx(wet)+u2hxn(wet)+g*h2x(wet)+g*h2xn(wet)+...
           +g*(hmean(wet)+hmeann(wet)).*bx(wet)+nu*(Del2u(wet)+Del2un(wet)));   
   hn(wet)=h(wet)-dt/2*(uhx(wet)+uhxn(wet));
   un(wet)=uhn(wet)./hn(wet);  
           
   hn=extrapolatingBC(hn+b,x,dry,wet,epsE,nn,dx)-b;
   un=extrapolatingBC(un,x,dry,wet,epsE,nn,dx);
       
          
  % Boundary conditions (Reflective)
   un(1)=-un(2); un(n)=-un(n-1);
   hn(1)=hn(2)+b(2)-b(1);  hn(n)=hn(n-1)+b(n-1)-b(n);
           
   h=hn; u=un;
      
   M(ii)=sum(h(wet)*dx);
   Herror(ii-1)=norm(h(wet)-hana(wet),2)./norm(hana(wet),2);
   Huerror(ii-1)=norm(u(wet).*h(wet)-uana(wet).*hana(wet),2);
   
end



figure
t1=linspace(0,t,length(M));
plot(t1,abs(M-M(1))/M(1),'k'); grid
fprintf('Max. mass error: %5.3g\n',max(abs(M-M(1))/M(1)));
xlabel('t'); ylabel('|M-M_0|/M_0'); set(gca,'FontSize',14);

figure
t2=linspace(0,t,length(Herror));
plot(t2,Herror,'k'); grid
fprintf('Max. h error: %5.3g\n',max(Herror));
xlabel('t'); ylabel('||h-h_{a}||_\infty/||h_{a}||_\infty'); set(gca,'FontSize',14);

figure
plot(t2,Huerror,'k'); grid
fprintf('Max. hu error: %5.3g\n',max(Huerror));
xlabel('t'); ylabel('||hu-(hu)_a||_\infty'); set(gca,'FontSize',14);

end



function P=extrapolatingBC(P,x,dry,wet,eps,nn,dx)

% P(dry)=rbfExtrapolationLocal(x(wet),P(wet),x(dry),eps,nn,dx);
P(dry)=rbfExtrapolationLocalPolynomial(x(wet),P(wet),x(dry),eps,nn,dx);
% P(dry)=rbfExtrapolationPolynomial(x(wet),P(wet),x(dry),1,dx);
   
end
