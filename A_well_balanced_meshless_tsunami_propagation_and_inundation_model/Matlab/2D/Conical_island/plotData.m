function F=plotData(ep)

m=200; n=200;
x=linspace(0,25,m);
y=linspace(0,28,n);

% use only every third point
%  x=x(1:3:end);m=length(x);
%  y=y(1:3:end);n=length(y);

dx=x(2)-x(1); dy=y(2)-y(1);
    t=20;
    dt=0.02;
[X,Y]=meshgrid(x,y);
steps=round(t/dt);
X=X'; Y=Y';

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
%         if(r<=2.2)
%             b(x==xc,y==yc)=0.625;            
%         elseif(r<=7.2)
%             b(x==xc,y==yc)=-r/8+0.9;
%         end                   
    end
end

d=0.32;
start=300;%2;%40;%100;% normally start at 2

delta=1/70;

% Case
%   ep=0.04 % Case A
% ep=0.09 % Case B
% ep=0.18; % Case C

load('testRunUp.mat')
load('gages.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show max runup
% figure
% 
%         surf(X,Y,b);
%         colormap copper
%         freezeColors;
%         hold on
%         colormap jet
%         surf(X,Y,hmax);
        
% compare to experimental data
Hmax=full(hmax);
% temp=Hmax;
 %  Hmax(Hmax-b<1/70)=nan;    
[theta,~]=cart2pol(X-12.96,Y-13.8);
B=Hmax(Hmax>delta)-b(Hmax>delta);
%    B=Hmax-b;
idx=find(abs(B)<0.01);
% idx1=find(abs(B(1:110,:))<=0.025);%0.01     %0.013% B<=0.000001
% idx2=find(abs(B(110:end,:))<=0.05);
T1=theta(Hmax>0);
H=Hmax(Hmax>0);
% T1=theta(1:110,:);
% T2=theta(110:end,:);
% H1=temp(1:110,:);
% H2=temp(110:end,:);

rotate=-1.5708;
load('measuredData.mat')

if(ep==0.04)
    measuredRunup=runupA;
    measuredGages=gagesA;
end
if(ep==0.09)
    measuredRunup=runupB;
    measuredGages=gagesB;
end
if(ep==0.18)
    measuredRunup=runupC;    
    measuredGages=gagesC;
end

% plot runup 
figure
polarplot(rad+rotate,110+((32+measuredRunup)-62.5)*(-4),'dblack');
rad2=linspace(0,2*pi,360);
hold on
polarplot(rad2,110+((100*d*ones(360))-62.5)*(-4),'black');
hold on
polarplot(T1(idx),110+(100*H(idx)-62.5)*(-4),'*black');
% polarplot(T1(idx1),110+(100*H1(idx1)-62.5)*(-4),'*black');
% hold on
% polarplot(T2(idx2),110+(100*H2(idx2)-62.5)*(-4),'*black');

set(gca,'fontsize',14);
% return
% plot gages
 

 figure
titl=[6,9,16,22];

l=[0,0,1,2];

for ii=1:4
%     figure
    subplot(1,4,ii)
    plot(time,measuredGages(:,ii),'--black')
    hold on
    if(ep==0.18)
        plot(linspace(20,t+20,steps),gages(:,ii)-d,'black')        
    end
    if(ep==0.04)
        plot(linspace(20,t+24,steps),gages(:,ii)-d,'black')
    end
    title(num2str(titl(ii)));
%     axis([22 42 -0.05 0.03]) % A
    axis([22 40 -0.045 0.1])
    grid
    set(gca,'fontsize',14);
    if(ii==1)
        xlabel('t')
        ylabel('h+b')
    end
end
