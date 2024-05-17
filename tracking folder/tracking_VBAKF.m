function res = tracking_VBAKF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
tao=3;N=10;
rou=1-exp(-4);
uapriv=(size(F_k,1)+1+tao);
Uapriv=tao*R_k;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1,uapriv,Uapriv]= VBAKF(x_k_1k_1,P_k_1k_1,uapriv,Uapriv,F_k,H_k,z_k,Q_k,R_k,N,tao,rou);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk,ukk,Ukk,D_Pk1k,D_R] = VBAKF(xkk,Pkk,ukk,Ukk,F,H,z,Q,R,N,tao1,rou)

%%%%%%
nz=size(z,1);

nx=size(xkk,1);

%%%%%Time update
xk1k=F*xkk;

Pk1k=F*Pkk*F'+Q;

tk1k=(nx+1+tao1);

Tk1k=tao1*Pk1k;

uk1k=rou*(ukk-nz-1)+nz+1;

Uk1k=rou*Ukk;

%%%%%Measurement update
xkk=xk1k;

Pkk=Pk1k;

for i=1:N
    
    %%%%%%%Update IW distribution parameters for Pk1k
    Ak=(xkk-xk1k)*(xkk-xk1k)'+Pkk;     %%%%
    
    tkk=tk1k+1;                        %%%%
    
    Tkk=Tk1k+Ak;                       %%%%
    
    E_i_Pk1k=(tkk-nx-1)*inv(Tkk);      %%%%
    
    %%%%%%%Update IW distribution parameters for Rk
    Bk=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';  %%%%
    
    ukk=uk1k+1;                        %%%%
    
    Ukk=Uk1k+Bk;                       %%%%
    
    E_i_Rk=(ukk-nz-1)*inv(Ukk);        %%%%
    
    %%%%%%%Update state estimate
    D_Pk1k=inv(E_i_Pk1k);
    
    D_R=inv(E_i_Rk);
    
    zk1k=H*xk1k;
    
    Pzzk1k=H*D_Pk1k*H'+D_R;
    
    Pxzk1k=D_Pk1k*H';
    
    Kk=Pxzk1k*inv(Pzzk1k);
    
    xkk=xk1k+Kk*(z-zk1k);
    
    Pkk=D_Pk1k-Kk*H*D_Pk1k;
end
end