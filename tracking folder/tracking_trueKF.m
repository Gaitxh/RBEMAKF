function res = tracking_trueKF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    Q_k = sim.Q(:,:,k);R_k = sim.R(:,:,k);
    [x_k_1k_1,P_k_1k_1]=KF(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk]=KF(xkk,Pkk,F,H,z,Q,R)

%%%%%Time update

xk1k=F*xkk;

Pk1k=F*Pkk*F'+Q;

%%%%%Measurement update

Pzzk1k=H*Pk1k*H'+R;

Pxzk1k=Pk1k*H';

Kk=Pxzk1k*inv(Pzzk1k);

xkk=xk1k+Kk*(z-H*xk1k);

Pkk=Pk1k-Kk*H*Pk1k;
end