function [Qk,Rk] = aaaQkEvaluate(xk_1k_1,Pk_1k_1,z_k,Fk,Hk)
x_kk_1 = Fk*xk_1k_1;
y__k = z_k-Hk*x_kk_1;
[n,m]=size(Hk);
H = Hk*Fk*Pk_1k_1*Fk'*Hk';
inv_H = inv(H);
wk = pinv(Hk'*Hk)*Hk'*y__k;
Qk = diag(wk.^2);
% [wk,err] = DLap_Estimate(Hk,y__k);
% Qk = diag(err);

rng('default')
Rk = abs(diag((y__k-Hk*wk-mvnrnd(zeros(1, n), H, 1)').^2));
rng('shuffle')

Pvk = Rk-Rk'/(H+Rk)*Rk;
phy_r= 1;
phy_q= 1;

count = 1;
while 1
    count = count  + 1;
    
    Omiga = inv_H-inv_H*Pvk*inv_H;
    Pwk = Qk-Qk*Hk'/(Hk*Qk*Hk'+inv(Omiga))*Hk*Qk;
    Ewk = Pwk*Hk'*Omiga*y__k;
    q(:,count) = (phy_q/4)*(sqrt(1+8*(Ewk.^2+diag(Pwk))/phy_q)-1);
%     phy_q= mean(q(:,count));
    Qk = diag(q(:,count));
    wk = Ewk;
    
    Pvk = Rk-Rk'/(H+Rk)*Rk;
    Evk = Pvk*inv_H*(y__k-Hk*wk);
    r(:,count) = (phy_r/4)*(sqrt(1+8*(Evk.^2+diag(Pvk))/phy_r)-1);
%     phy_r= mean(r(:,count));
    Rk = diag(r(:,count));
    vk = Evk;
    
    
    target(count) = .5*(norm(q(:,count)-q(:,count-1),1)/norm(q(:,count-1),1)...
        +norm(r(:,count)-r(:,count-1),1)/norm(r(:,count-1),1));
    if target(count) < 1e-3 || count > 30
        break;
    end
end
end