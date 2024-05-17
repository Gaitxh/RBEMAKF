function res = tracking_KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    N = 10;v = 5;tao_P = 5;
    [x_k_1k_1,P_k_1k_1]=stu2_kf(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k,N,v,v,tao_P);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk]=stu2_kf(xkk,Pkk,F,H,z,Q,R,N,v1,v2,tao)

%%%%%%Set up
nz=size(z,1);

nx=size(xkk,1);

%%%%%Time update
xk1k=F*xkk;

Pk1k=F*Pkk*F'+Q;

%%%%%Measurement update
xkk=xk1k;

Pkk=Pk1k;

t0=(nx+1+tao);

T0=tao*Pk1k;

E_i_Pk1k=(t0-nx-1)*inv(T0);

for i=1:N
    
    %%%%%%%Update the distributions of \xi and \lambda
    Dk1=(xkk-xk1k)*(xkk-xk1k)'+Pkk;
    
    Dk2=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';
    
    gama1=trace(Dk1*E_i_Pk1k);
    
    gama2=trace(Dk2*inv(R));
    
    E_kasai=(v1+nx)/(v1+gama1);
    
    E_lamda=(v2+nz)/(v2+gama2);
    
    %%%%%%%Update the distribution of Pk1k
    tk1k=t0+1;
    
    Tk1k=T0+E_kasai*Dk1;
    
    E_i_Pk1k=(tk1k-nx-1)*inv(Tk1k);
    
    %%%%%%%Update the distribution of state vector
    D_Pk1k=inv(E_i_Pk1k)/E_kasai;
    
    D_R=R/E_lamda;
    
    zk1k=H*xk1k;
    
    Pzzk1k=H*D_Pk1k*H'+D_R;
    
    Pxzk1k=D_Pk1k*H';
    
    Kk=Pxzk1k*inv(Pzzk1k);
    
    xkk=xk1k+Kk*(z-zk1k);
    
    Pkk=D_Pk1k-Kk*H*D_Pk1k;
    
end
end