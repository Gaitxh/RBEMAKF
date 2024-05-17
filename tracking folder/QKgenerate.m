function [Qk,Rk] = QKgenerate(Hk,Fk,Pk_1k_1,zk,xk1k1)
xk1k = Fk*xk1k1;
y__k = zk-Hk*xk1k;
[n,m] = size(Hk);
H = Hk*Fk*Pk_1k_1*Fk'*Hk';
inv_H = inv(H);
wk = pinv(Hk'*Hk)*Hk'*y__k;

HkFkxk = utchol(H)*randn(n,1);
res = y__k-Hk*wk-HkFkxk;
Qk = diag(wk.^2);
Rk = diag(res.^2);
end