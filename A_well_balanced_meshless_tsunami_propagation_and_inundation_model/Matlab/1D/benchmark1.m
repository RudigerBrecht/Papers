function benchmark1

rng(1);
close all

   
t=20;
dt=0.0015;

steps=round(t/dt);

g=9.81;

n=50;
x=linspace(0,1,n)';
dx=x(2)-x(1);

% Threshold for dry areas
delta=0.0025;

% Bottom topography
bs=zeros(n,1);
rm=0.4;
a=1.2;

for i=1:n
    r=abs(x(i)-0.5);
   if(r<0.4)
       bs(i)=a*(exp(-0.5/(rm^2-r^2)))/(exp(-0.5/rm^2));
   end
end

bn=zeros(n,1);
a=[0.1,0.2,0.3];
p=[1.6,3.2,0.5]';
for i=1:n
    bn(i)=bs(i)+sum(a*sin(16*[1;2;3]*pi*x(i)+p));
end

b=bs;

% Initial conditions
u=zeros(n,1);


h=ones(n,1);
% h(n/4-2:n/4+2)=h(n/4-2:n/4+2)+0.2;
% for ii=1:3
%   h(ix)=(h(ix+1)+h(ix)+h(ix-1))/3;
% end
h=h-b;
h0=h;



% Differentation matrix
eps=0.8;%2; 
nn=3;
[Dx,Del2]=rbfFD1GApoly(x',eps,nn,dx);


% Conserved quantities
M=zeros(steps,1);
M(1)=sum(h*dx);

Herror=zeros(steps,1);

% Time stepping
un=u; hn=h; uhn=u.*h;

nu=-0.01;epsE=1;
for ii=2:steps
    
   fprintf('Computing %ith step out of %i steps.\n',ii,steps);
   

   wet=h>delta; dry=h<=delta;  
   
   % Derivatives for conservative form
   u2hx=Dx*(u.^2.*h); % h2x=0.5*Dx*(h.^2);
   uhx=Dx*(h.*u); bx=Dx*b; %hx=Dx*h;
   Del2u=Del2*(h.*u);
      
   % Averaging matrix
%    h2x=computeDxh2(MDx,h); 
           
   % Euler forward step
   Ma=Dx; Ma(Dx~=0)=1./nn;
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
    Ma=Dx; Ma(Dx~=0)=1./nn;
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
   Herror(ii-1)=norm(h(wet)-h0(wet),inf)./norm(h0(wet),inf);
   Huerror(ii-1)=norm(u(wet).*h(wet),inf);
   
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
