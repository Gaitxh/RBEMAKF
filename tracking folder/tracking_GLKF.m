function res = tracking_KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
Lambda = 0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    if sim.flag == 1
        [Q_k,R_k] = QKgenerate(H_k,F_k,P_k_1k_1,z_k,x_k_1k_1);
    elseif sim.flag == 2
        Q_k = sim.Q(:,:,k);R_k = sim.R(:,:,k);
    end
    [x_k_1k_1,P_k_1k_1,R_k]=GLKF(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k,Lambda,sim);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [x_k_1k_1,P_k_1k_1,R_k_1]=GLKF(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k_1,Lambda,sim)

    xkk_1 = F_k*x_k_1k_1;
    E_k = z_k-H_k*xkk_1;
    Pkk_1 = F_k*P_k_1k_1*F_k'+Q_k;
    R_k = (1-Lambda)*R_k_1 + eye(sim.n)*Lambda*(E_k'*E_k)/(sim.length-1);
    R_k_1 = R_k;
    y_k = z_k-H_k*xkk_1;
    S_k = H_k*Pkk_1*H_k'+trace(R_k)*eye(sim.n);
    K_k = Pkk_1*H_k'/S_k;
    x_k_1k_1 = x_k_1k_1 + K_k*y_k;
    P_k_1k_1 = Pkk_1-K_k*H_k*Pkk_1;
end