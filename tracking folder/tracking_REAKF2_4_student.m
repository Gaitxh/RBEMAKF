function res = tracking_EAKF2(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Qk = sim.Q0;
Rk = sim.R0;
vq = 5;vr = 5;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1,Qk,Rk]=adaptivefc(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Qk,Rk,vq,vr);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk,Pwk,Pvk] = ...
    adaptivefc(zk,Hk,Fk,xk1k1,Pk1k1,Qk,Rk,vq,vr)
xk1k = Fk*xk1k1;
y__k = zk-Hk*xk1k;

[n,m] = size(Hk);
H = Hk*Fk*Pk1k1*Fk'*Hk';
inv_H = inv(H);
Pvk = Rk-Rk'/(H+Rk)*Rk;
q = diag(Qk);
r = diag(Rk);
for i = 2:30
    Omiga = inv_H-inv_H*Pvk*inv_H;
    Pwk = Qk-Qk*Hk'/(Hk*Qk*Hk'+inv(Omiga))*Hk*Qk;
    Ewk = Pwk*Hk'*Omiga*y__k;
    q = (Ewk.^2+diag(Pwk)+vq-2)./(vq+3);
    Qk = diag(q);    
    wk = Ewk;
    
    Pvk = Rk-Rk'/(H+Rk)*Rk;
    Evk = Pvk*inv_H*(y__k-Hk*wk);
    r = (Evk.^2+diag(Pvk)+vr-2)./(vr+3);
    Rk = diag(r);
    vk = Evk;
    
    Pk1k=Fk*Pk1k1*Fk'+Pwk;
    Pzzk1k=Hk*Pk1k*Hk'+Pvk;
    Pxzk1k=Pk1k*Hk';
    Kk=Pxzk1k*inv(Pzzk1k);
    xkk(:,i)=xk1k+Kk*(zk-Hk*xk1k);
    Pkk=Pk1k-Kk*Hk*Pk1k;
    target_fuc(i)=norm(xkk(:,i)-xkk(:,i-1),1)/norm(xkk(:,i-1),1);
    if target_fuc(i)<1e-8
        break;
    end
end
xkk = xkk(:,end);
end
