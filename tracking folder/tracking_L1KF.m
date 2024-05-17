function res = tracking_L1KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1]=l1kf(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk,Pk1k]=l1kf(xk1k1,Pk1k1,A,H,z,Q,R)

Lambda = pinv(H'*H)*H'*z-A*xk1k1;
Ck = Lambda'*sign(Lambda);
xk1k=A*xk1k1;
Pk1k=A*Pk1k1*A'+Ck'*Q*Ck;
Kk=(Pk1k*H')*inv(H*Pk1k*H'+R);
xkk=xk1k+Kk*(z-H*xk1k);
Pkk=Pk1k-Kk*H*Pk1k;
end