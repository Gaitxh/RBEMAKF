function res = tracking_BayesOutputKF(sim)

H_k = sim.H;
F_k = sim.F;
m = size(F_k,1); n = size(H_k,1);
% inivalue defined
%%set all of them to zero vectoes and matrices of appropriate dimension
zk_1k_1 = sim.inix;
Pzkk_1 = sim.iniP;
vk_1 = 0;rouk_1 = 0;
% consider as the identity matrices of appropriate dimensions multiplied by small scalars
Omigak_1 = eye(size(F_k,1))*1e-1;
Sigmak_1 = eye(size(H_k,1))*1e-1;
Rk_1= Sigmak_1/(vk_1+n+1);

zkk_1 = F_k*zk_1k_1;

for k = 1 : sim.length
    z_k = sim.z(:,k);
    dk=z_k;G = H_k;J=H_k;A=F_k;B=F_k;N0=n;N=m;
    [Rk_1,Pzkk_1,zkk_1,Sigmak_1,vk_1,Omigak_1,rouk_1,zk_1k_1] = ...
    adaptive1_KF(Rk_1,dk,G,Pzkk_1,zkk_1,Sigmak_1,vk_1,N0,Omigak_1,rouk_1,A,zk_1k_1,N);
    rec_state(:,k) = zk_1k_1;
    rec_data(:,k) = H_k*zk_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [Rk,Pzk1k,zk1k,Sigmak,vk,Omigak,rouk,zkk] = ...
    adaptive1_KF(Rk_1,dk,G,Pzkk_1,zkk_1,Sigmak_1,vk_1,N0,Omigak_1,rouk_1,A,zk_1k_1,N)
%%state estimation
Gzk = Pzkk_1*G'/(G*Pzkk_1*G'+Rk_1);
zkk = zkk_1+Gzk*(dk-G*zkk_1);
Pzkk = Pzkk_1-Gzk*G*Pzkk_1;
%%update the covariance matrix of observation noise
Sigmak=Sigmak_1+(dk-G*zkk)*(dk-G*zkk)';
vk = vk_1 + 1;
Rk=Sigmak/(vk+N0+1);
%%update the covariance matrix of process noise
Omigak = Omigak_1 + (zkk-A*zk_1k_1)*(zkk-A*zk_1k_1)';
rouk = rouk_1 + 1;
Qk = Omigak/(rouk+2*N+1);
%%predict the state
zk1k = A*zkk;
Pzk1k = A*Pzkk*A'+Qk;
end