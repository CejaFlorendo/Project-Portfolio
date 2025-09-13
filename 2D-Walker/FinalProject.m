%2D walker 
clear;close all;

t_sim=10;

dt=0.02;

%leg geometry
L1=0.4;      %thigh [m]

L2=0.3;      %shank  [m]

H_hip=0.69;  %fixed hip height [m]

%masses/physics
M_hip=20;    %hip mass [kg]

M_thigh=1.0; %thigh mass [kg]

M_shank=1.0; %shin mass [kg]

g=9.81;

%joint PD values
Kp=50;

Kd=5;

%ground contact (stance leg only)
k_ground=40000;  %vertical stiffness [N/m]

u=0.8;           %CoF

%extra damping at stance to prevent slipping
c_t=80;          %[N·s/m] along x

%gait timing 
stride_freq=1;

omega=2*pi*stride_freq;

%desired joint trajectories
hip_swing=pi/9;

knee_swing=pi/6;

S=@(t) sin(omega*t);

C=@(t) cos(omega*t);

%hip targets
hipAngDesL=@(t) -hip_swing*S(t);

hipAngVelDesL=@(t) -hip_swing*omega*C(t);

hipAngDesR=@(t) -hip_swing*S(t+pi/omega);  

hipAngVelDesR=@(t) -hip_swing*omega*C(t+pi/omega);

%knee targets
kneeAngDesL=@(t) (S(t)>=0).*(knee_swing.*sin(pi*S(t)));                          %right stance left swings

kneeAngVelDesL=@(t) (S(t)>=0).*(knee_swing.*cos(pi*S(t)).*(pi*omega.*C(t)));

kneeAngDesR=@(t) (S(t)<0).*(knee_swing.*sin(pi*(-S(t))));                        %left stance right swings

kneeAngVelDesR=@(t) (S(t)<0).*(knee_swing.*cos(pi*(-S(t))).*(pi*(-omega).*C(t)));

%target hip speed
v_target=0.3;   %+right, −left [m/s]

Kp_v=120;    %hip vel P gain

Kd_v=0;

%initial state
X0=zeros(10,1);

X0(1)=0.0;                 %hip_x inital    

X0(2)=0.0;                 %hip_v initial

X0(3)=hipAngDesL(0);

X0(4)=kneeAngDesL(0);

X0(5)=hipAngVelDesL(0); 

X0(6)=kneeAngVelDesL(0);

X0(7)=hipAngDesR(0); 

X0(8)=kneeAngDesR(0);

X0(9)=hipAngVelDesR(0); 

X0(10)=kneeAngVelDesR(0);

%params
par=struct('L1',L1,'L2',L2,'hip_h',H_hip,'m_hip',M_hip,'m_th',M_thigh,'m_sh',M_shank,'g',g,...
           'Kp',Kp,'Kd',Kd,'k_ground',k_ground,'mu',u,'omega',omega,'hipAngDesL',hipAngDesL,'hipAngVelDesL',hipAngVelDesL,'hipAngDesR',hipAngDesR,'hipAngVelDesR',hipAngVelDesR,...
           'kneeAngDesL',kneeAngDesL,'kneeAngVelDesL',kneeAngVelDesL,'kneeAngDesR',kneeAngDesR,'kneeAngVelDesR',kneeAngVelDesR,'Kp_v',Kp_v,'Kd_v',Kd_v,'v_des',v_target,'c_t',c_t);

%%integrate
opts=odeset('RelTol',1e-6,'AbsTol',1e-8);

[t_out,X]=ode45(@(tt,xx) f_ode(tt,xx,par),[0 t_sim],X0,opts);

hip_x=X(:,1);

thL=X(:,3); 

delL=X(:,4);

thR=X(:,7); 

delR=X(:,8);

%%animation
fig=figure('Color',[0 0 0],'InvertHardcopy','off');
ax=axes('Parent',fig,'Color',[0.06 0.06 0.06],'XColor','w','YColor','w');
hold(ax,'on'); axis(ax,'equal');
xlabel(ax,'X (m)','Color','w');
ylabel(ax,'Y (m)','Color','w');
title(ax,'2D walker','Color','w');

cam_half=2.0;      

ymin=-0.2; 

ymax=1.0;

step=max(1,round(dt/mean(diff(t_out))));

for k=1:step:length(t_out)

    cla(ax);

    hx=hip_x(k); hy=H_hip;

    %follow cam
    set(ax,'XLim',[hx-cam_half,hx+cam_half],'YLim',[ymin ymax]);

    %ground across current limits
    xl=get(ax,'XLim');

    plot(ax,xl,[0 0],'w-','LineWidth',2); hold(ax,'on');

    %left leg (blue)
    kxL=hx+L1*sin(thL(k)); 

    kyL=hy-L1*cos(thL(k));

    fxL=kxL+L2*sin(thL(k)+delL(k));

    fyL=kyL-L2*cos(thL(k)+delL(k));

    %right leg (red)
    kxR=hx+L1*sin(thR(k));

    kyR=hy-L1*cos(thR(k));

    fxR=kxR+L2*sin(thR(k)+delR(k));

    fyR=kyR-L2*cos(thR(k)+delR(k));

    %clamp stance foot to y>=0
    s=sin(omega*t_out(k));

    if s<0,fyL=max(0,fyL);

    else,

        fyR=max(0,fyR);

    end

    plot(ax,[hx kxL fxL],[hy kyL fyL],'-o','Color',[0 0.5 1],...
         'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[0 0.5 1]);
    plot(ax,[hx kxR fxR],[hy kyR fyR],'-o','Color',[1 0.3 0.3],...
         'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',[1 0.3 0.3]);

    drawnow;
end

%plot 
figure('Color',[0 0 0],'InvertHardcopy','off');
ax2=axes('Color',[0.06 0.06 0.06],'XColor','w','YColor','w');
plot(ax2,t_out,hip_x,'w','LineWidth',1.5); grid(ax2,'on');
xlabel(ax2,'t (s)','Color','w'); ylabel(ax2,'hip x (m)','Color','w');
title(ax2,'Hip translation','Color','w');

vxL=zeros(size(t_out));
vyL=vxL;

vxR=zeros(size(t_out)); 
vyR=vxR;

for k=1:length(t_out)
    hx=X(k,1); 

    hv=X(k,2);

    thl=X(k,3); 

    del_l=X(k,4); 

    dthl=X(k,5); 

    ddel_l=X(k,6);

    thr=X(k,7); 

    del_r=X(k,8); 

    dthr=X(k,9); 

    ddel_r=X(k,10);

    [~,~,xdL,ydL,~,~]=foot_kin(hx,H_hip,thl,del_l,hv,dthl,ddel_l,L1,L2);

    [~,~,xdR,ydR,~,~]=foot_kin(hx,H_hip,thr,del_r,hv,dthr,ddel_r,L1,L2);

    vxL(k)=xdL; vyL(k)=ydL;

    vxR(k)=xdR; vyR(k)=ydR;
end

%foot velocity plots
figV=figure('Color',[0 0 0],'InvertHardcopy','off'); 
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

axv1=nexttile; set(axv1,'Color',[0.06 0.06 0.06],'XColor','w','YColor','w'); hold(axv1,'on');
plot(axv1,t_out,vxL,'Color',[0 0.5 1],'LineWidth',1.5);
plot(axv1,t_out,vxR,'Color',[1 0.3 0.3],'LineWidth',1.5);
grid(axv1,'on'); ylabel(axv1,'foot v_x (m/s)','Color','w');
title(axv1,'Foot horizontal velocity','Color','w'); legend(axv1,{'Left','Right'},'TextColor','w','Location','best');

axv2=nexttile; set(axv2,'Color',[0.06 0.06 0.06],'XColor','w','YColor','w'); hold(axv2,'on');
plot(axv2,t_out,vyL,'Color',[0 0.5 1],'LineWidth',1.5);
plot(axv2,t_out,vyR,'Color',[1 0.3 0.3],'LineWidth',1.5);
grid(axv2,'on'); xlabel(axv2,'t (s)','Color','w'); ylabel(axv2,'foot v_y (m/s)','Color','w');
title(axv2,'Foot vertical velocity','Color','w'); legend(axv2,{'Left','Right'},'TextColor','w','Location','best');

%knee angle plots
figK=figure('Color',[0 0 0],'InvertHardcopy','off');
axk=axes('Color',[0.06 0.06 0.06],'XColor','w','YColor','w'); hold(axk,'on');
plot(axk,t_out,rad2deg(delL),'Color',[0 0.5 1],'LineWidth',1.5);
plot(axk,t_out,rad2deg(delR),'Color',[1 0.3 0.3],'LineWidth',1.5);
grid(axk,'on');
xlabel(axk,'t (s)','Color','w'); ylabel(axk,'knee angle (deg)','Color','w');
title(axk,'Knee flexion angle','Color','w'); legend(axk,{'Left','Right'},'TextColor','w','Location','best');

function dX=f_ode(t,X,par)

    hip_x=X(1); hip_v=X(2);

    thL=X(3); delL=X(4); dthL=X(5); ddelL=X(6);

    thR=X(7); delR=X(8); dthR=X(9); ddelR=X(10);

    s=sin(par.omega*t);

    L_stance=(s<0);

    
    thL_d=par.hipAngDesL(t); dthL_d=par.hipAngVelDesL(t);

    thR_d=par.hipAngDesR(t); dthR_d=par.hipAngVelDesR(t);

    delL_d=par.kneeAngDesL(t); ddelL_d=par.kneeAngVelDesL(t);

    delR_d=par.kneeAngDesR(t); ddelR_d=par.kneeAngVelDesR(t);

    [ML,CL,GL]=leg_dyn_mats(thL,delL,dthL,ddelL,par.m_th,par.m_sh,par.L1,par.L2,par.g);

    [MR,CR,GR]=leg_dyn_mats(thR,delR,dthR,ddelR,par.m_th,par.m_sh,par.L1,par.L2,par.g);

    tauL_pd=[par.Kp*(thL_d-thL)+par.Kd*(dthL_d-dthL); par.Kp*(delL_d-delL)+par.Kd*(ddelL_d-ddelL)];

    tauR_pd=[par.Kp*(thR_d-thR)+par.Kd*(dthR_d-dthR); par.Kp*(delR_d-delR)+par.Kd*(ddelR_d-ddelR)];

    [~,yL,xdL,~,JxL,JyL]=foot_kin(hip_x,par.hip_h,thL,delL,hip_v,dthL,ddelL,par.L1,par.L2);

    [~,yR,xdR,~,JxR,JyR]=foot_kin(hip_x,par.hip_h,thR,delR,hip_v,dthR,ddelR,par.L1,par.L2);

    FyL=0; FyR=0; FxL=0; FxR=0;

    %vertical contact
    if L_stance

        if yL<0, FyL=par.k_ground*(-yL);

        end

        FxL=par.Kp_v*(par.v_des-hip_v);
        
        FxL=FxL-par.c_t*xdL;
        
        if FyL>0, FxL=max(-par.mu*FyL,min(par.mu*FyL,FxL)); 

        else, 

            FxL=0; 

        end

    else
        
        if yR<0, FyR=par.k_ground*(-yR); 
        
        end
        
        FxR=par.Kp_v*(par.v_des-hip_v);
        
        FxR=FxR-par.c_t*xdR;
        
        if FyR>0, FxR=max(-par.mu*FyR,min(par.mu*FyR,FxR)); 
        
        else, FxR=0; 
        
        end
    end

    tauL=tauL_pd+[JxL(1)*FxL+JyL(1)*FyL; JxL(2)*FxL+JyL(2)*FyL];

    tauR=tauR_pd+[JxR(1)*FxR+JyR(1)*FyR; JxR(2)*FxR+JyR(2)*FyR];

    ddqL=ML\(tauL-CL*[dthL; dthL+ddelL]-GL);

    ddqR=MR\(tauR-CR*[dthR; dthR+ddelR]-GR);

    hip_a=(FxL+FxR)/par.m_hip;

    dX=zeros(10,1);

    dX(1)=hip_v;

    dX(2)=hip_a;

    dX(3)=dthL;

    dX(4)=ddelL;

    dX(5)=ddqL(1);

    dX(6)=ddqL(2);

    dX(7)=dthR;

    dX(8)=ddelR;

    dX(9)=ddqR(1);

    dX(10)=ddqR(2);
end

function [M,C,G]=leg_dyn_mats(theta1,delta,dtheta1,ddelta,m_th,m_sh,L1,L2,g)
    s=sin(delta); 
    
    c=cos(delta);

    M=[(m_th+m_sh)*L1^2+m_sh*L2^2+2*m_sh*L1*L2*c, m_sh*L2^2+m_sh*L1*L2*c; m_sh*L2^2+m_sh*L1*L2*c, m_sh*L2^2];

    dtheta2=dtheta1+ddelta;

    C=[-m_sh*L1*L2*s*dtheta2,-m_sh*L1*L2*s*(dtheta1+dtheta2); m_sh*L1*L2*s*dtheta1,0];

    theta2=theta1+delta;

    G=[(m_th+m_sh)*g*L1*sin(theta1)+m_sh*g*L2*sin(theta2); m_sh*g*L2*sin(theta2)];
end

function [x,y,xd,yd,Jx,Jy]=foot_kin(hip_x,hip_h,theta1,delta,hip_v,dtheta1,ddelta,L1,L2)

    th2=theta1+delta;

    x=hip_x+L1*sin(theta1)+L2*sin(th2);

    y=hip_h-L1*cos(theta1)-L2*cos(th2);

    xd=hip_v+L1*cos(theta1)*dtheta1+L2*cos(th2)*(dtheta1+ddelta);

    yd=L1*sin(theta1)*dtheta1+L2*sin(th2)*(dtheta1+ddelta);

    Jx=[L1*cos(theta1)+L2*cos(th2); L2*cos(th2)];

    Jy=[L1*sin(theta1)+L2*sin(th2); L2*sin(th2)];
end
