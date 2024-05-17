function res = tracking_MLKF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
phy_k = 1;
rou = size(H_k,1)+3;
N = 10;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1] = MLKF(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Q_k,R_k,size(F_k,1),size(H_k,1),rou,phy_k,N);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [x_k_1k_1,P_k_1k_1] = MLKF(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Q,R,m,n,rou,phy_k,N)
%%time update
xkk_1 = F_k*x_k_1k_1;
Pkk_1 = F_k*P_k_1k_1*F_k'+Q;
%%Inilization
xkk(:,1)=xkk_1;Pkk(:,:,1)=Pkk_1;
V=rou*R;
EinvRk(:,:,1)=(rou+1)*eye(n)/V;
Einvpai(1) = 1;
ak = 1-.5*n; ck = phy_k;
for i = 1:N
    B_k(:,:,i+1)=(z_k-H_k*xkk(:,i))*(z_k-H_k*xkk(:,i))'+H_k*Pkk(:,:,i)*H_k';
    
    bk(i+1)=.5*trace(B_k(:,:,i+1)*EinvRk(:,:,i));
    Einvpai(i+1) = sqrt(ck/bk(i+1))*besselk(ak-1,2*sqrt(bk(i+1)*ck))/besselk(ak,2*sqrt(bk(i+1)*ck));
    if isnan(Einvpai(i+1)) == 1 || isinf(Einvpai(i+1)) == 1; Einvpai(i+1) = 1;end
    
    rou(i+1)=rou(1)+1; V(:,:,i+1) = Einvpai(i+1)*B_k(:,:,i+1) + V(:,:,1);
    EinvRk(:,:,i+1)=(rou(i+1))*eye(n)/V(:,:,i+1);
    
    R(:,:,i+1) = inv(EinvRk(:,:,i+1));
    Pkk_1(:,:,i+1) = Pkk_1(:,:,1);
    KalmanGain(:,:,i+1) = Pkk_1(:,:,i+1)*H_k'/(H_k*Pkk_1(:,:,i+1)*H_k'+R(:,:,i+1));
    xkk(:,i+1) = xkk_1+KalmanGain(:,:,i+1)*(z_k-H_k*xkk_1);
    Pkk(:,:,i+1) = (eye(m)-KalmanGain(:,:,i+1)*H_k)*Pkk_1(:,:,i+1);
    Target(i) = norm(xkk(:,i+1)-xkk(:,i),1)/norm(xkk(:,i),1);
    if Target(i)<1e-4 || i == N
        x_k_1k_1 = xkk(:,i+1);
        P_k_1k_1 = Pkk(:,:,i+1);
        break;
    end
end
end