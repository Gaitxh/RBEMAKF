function res = tracking_KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    N=50;              %%%%%%最大变分迭代次数
    sigma=5;           %%%%%%高斯核带宽
    v_sq=5;            %%%%%%自由度参数
    tao_P=2;
    tao_R=2;
    epsilon=1e-8;
    [x_k_1k_1,P_k_1k_1]=aorkf(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k,sigma,v_sq,tao_P,tao_R,N,epsilon,3);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk]=aorkf(xkk,Pkk,F,H,z,Q,R,sigma,v,tao_P,tao_R,N,epsilon,flag)

%%%%%%准备
nz=size(z,1);

nx=size(xkk,1);

%%%%%时间更新
xk1k=F*xkk;

Pk1k=F*Pkk*F'+Q;

%%%%%量测更新
xkk=xk1k;

Pkk=Pk1k;

lamda1=1;

lamda2=1;

%%%%%
uk1k=tao_R;

Uk1k=tao_R*R;

tk1k=tao_P;

Tk1k=tao_P*Pk1k;

E_Pk1k=Pk1k;

E_R=R;

for i=1:N
    
    %%%%%%%
    xkk_i=xkk;
    
    %%%%%%%计算状态估计
    D_Pk1k=E_Pk1k/lamda1;
    
    D_R=E_R/lamda2;
    
    zk1k=H*xk1k;
    
    Pzzk1k=H*D_Pk1k*H'+D_R;
    
    Pxzk1k=D_Pk1k*H';
    
    Kk=Pxzk1k*inv(Pzzk1k);
    
    xkk=xk1k+Kk*(z-zk1k);
    
    Pkk=D_Pk1k-Kk*H*D_Pk1k;
    
    %%%%%%%%判定
    td=norm(xkk-xkk_i)/norm(xkk);
    
    if td<=epsilon
        break;
    end
    
    %%%%%%%计算辅助参数
    Dk1=(xkk-xk1k)*(xkk-xk1k)'+Pkk;
    
    Dk2=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';
    
    %%%%%%%估计E_Pk1k
    tkk=tk1k+0.5;
    Tkk=Tk1k+0.5*lamda1*Dk1;
    E_Pk1k=Tkk/tkk;
    
    %%%%%%%估计E_R
    ukk=uk1k+0.5;
    Ukk=Uk1k+0.5*lamda2*Dk2;
    E_R=Ukk/ukk;
    
    %%%%%%%计算辅助参数
    gama1=trace(Dk1*inv(E_Pk1k));
    gama2=trace(Dk2*inv(E_R));
    
    %%%%%%%
    [lamda1,lamda2]=fun_ap(gama1,gama2,sigma,v,flag,nx,nz);
    
end
end