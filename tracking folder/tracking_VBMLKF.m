function res = tracking_VBMLKF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
tao = 5; rou = 1-exp(-4);
u_k_1k_1 = sim.n+3;U_k_1k_1 = u_k_1k_1*R_k;
phy_k = 1;N = 10;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    if sim.flag == 1
        [Q_k,R_k] = QKgenerate(H_k,F_k,P_k_1k_1,z_k,x_k_1k_1);
    end
    [x_k_1k_1,P_k_1k_1,u_k_1k_1,U_k_1k_1] = ...
        VBLKF(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Q_k,sim.m,sim.n,tao,rou,u_k_1k_1,U_k_1k_1,phy_k,N);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [x_k_1k_1,P_k_1k_1,u_k_1k_1,U_k_1k_1] = ...
    VBLKF(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Q_k,m,n,tao,rou,u_k_1k_1,U_k_1k_1,phy_k,N)
% Time update
x_kk_1 = F_k*x_k_1k_1;
P_kk_1 = F_k*P_k_1k_1*F_k'+Q_k;
% Inilization
x_kk(:,1) = x_kk_1; P_kk(:,:,1) = P_kk_1;
t_kk_1 = m+tao+1; T_kk_1 = tao*P_kk_1;
u_kk_1 = rou*(u_k_1k_1-n-1)+n+1; U_kk_1 = rou*U_k_1k_1;
E_invPkk_1(:,:,1)=(t_kk_1-m-1)*eye(m)*inv(T_kk_1);E_invRk(:,:,1)=(u_kk_1-n-1)*eye(n)*inv(U_kk_1);
E_invpai(1) = 1;
a_k = 1-.5*n; c_k = phy_k;
for i = 1:N
    A_k(:,:,i+1)=(z_k-H_k*x_kk(:,i))*(z_k-H_k*x_kk(:,i))'+H_k*P_kk(:,:,i)*H_k';
    B_k(:,:,i+1)=P_kk(:,:,i)+(x_kk(:,i)-x_kk_1)*(x_kk(:,i)-x_kk_1)';
    
    b_k(i+1)=.5*trace(A_k(:,:,i+1)*E_invRk(:,:,i));
    E_invpai(i+1) = sqrt(c_k/b_k(i+1))*besselk(a_k-1,2*sqrt(b_k(i+1)*c_k))/besselk(a_k,2*sqrt(b_k(i+1)*c_k));
    if isnan(E_invpai(i+1)) == 1 || isinf(E_invpai(i+1)) == 1; E_invpai(i+1) = 1;end
    
    u_k(i+1)=u_kk_1+1; U_k(:,:,i+1) = U_kk_1+E_invpai(i+1)*A_k(:,:,i+1);
    E_invRk(:,:,i+1)=(u_k(i+1)-n-1)*eye(n)*inv(U_k(:,:,i+1));
    
    t_k(i+1)=t_kk_1+1; T_k(:,:,i+1) = T_kk_1+B_k(:,:,i+1);
    E_invPkk_1(:,:,i+1)=(t_k(i+1)-m-1)*eye(m)*inv(T_k(:,:,i+1));
    
    R_k(:,:,i+1) = (E_invpai(i+1)\rou*U_k_1k_1+A_k(:,:,i+1))/(rou*(u_k_1k_1-n-1)+1);
    P_kk_1(:,:,i+1) = (tao*P_kk_1(:,:,1)+B_k(:,:,i+1))/(tao+1);
    KalmanGain(:,:,i+1) = P_kk_1(:,:,i+1)*H_k'*inv(H_k*P_kk_1(:,:,i+1)*H_k'+R_k(:,:,i+1));
    x_kk(:,i+1) = x_kk_1+KalmanGain(:,:,i+1)*(z_k-H_k*x_kk_1);
    P_kk(:,:,i+1) = (eye(m)-KalmanGain(:,:,i+1)*H_k)*P_kk_1(:,:,i+1);
    Target(i) = norm(x_kk(:,i+1)-x_kk(:,i),1)/norm(x_kk(:,i),1);
    if Target(i)<1e-4 || i == N
        x_k_1k_1 = x_kk(:,i+1);
        P_k_1k_1 = P_kk(:,:,i+1);
        u_k_1k_1 = u_k(i+1);
        U_k_1k_1 = U_k(:,:,i+1);
        break;
    end
end
end