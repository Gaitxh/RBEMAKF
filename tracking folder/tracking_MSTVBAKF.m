function res = tracking_MSTVBAKF(sim)

H = sim.H;
xmstvbakf = sim.inix;
Pmstvbakf = sim.iniP;
F = sim.F;
Q0 = sim.Q0;
R0 = sim.R0;
Rmstvbakf = R0;
N=10;
rou=1-exp(-4);
nz=sim.n;
tao_R=3;
umstvbakf=(nz+1+tao_R);
Umstvbakf=tao_R*R0;
Sk2=0;
for k = 1 : sim.length
    z = sim.z(:,k);
    
    [xmstvbakf,Pmstvbakf,umstvbakf,Umstvbakf,~,Rmstvbakf,Sk2]=mstvbakf(xmstvbakf,Pmstvbakf,umstvbakf,Umstvbakf,F,H,z,Q0,Rmstvbakf,N,rou,k,Sk2,sim.n,sim.m);
    rec_state(:,k) = xmstvbakf;
    rec_data(:,k) = H*xmstvbakf;
end
res.x = rec_state;
res.z = rec_data;
end

function [xkk,Pkk,ukk,Ukk,Pk1k,D_R,Sk]=mstvbakf(xkk,Pkk,ukk,Ukk,F,H,z,Q,R,N,rou,i,Sk,n,m)
%
nz=size(z,1);
yip=0.85;
beta=0.4;
% Time update
xk1k=F*xkk;
v = z - H * xk1k;
dk=(1-yip)/(1-yip^i);
Sk=dk*v*v'+(1-dk)*Sk;
Mk=H*F*Pkk*F'*H';
Nk=Sk-H*Q*H'-beta*R;
trace_Mk=0;
for ii=1:m
    if ii <= m/2
        a(ii,1) = 1;
    elseif ii > m/2
        a(ii,1) = 1.2;
    end
end
% a(1,1)=1;a(2,1)=1;a(3,1)=1.2;a(4,1)=1.2;
for ii=1:n
    trace_Mk=trace_Mk+a(ii,1)*Mk(ii,ii);
end
lemda=max(1,trace(Nk)/trace_Mk);
for ii=1:m
    Lemda(ii,ii)=a(ii,1)*lemda;
end
Pk1k=Lemda*(F*Pkk*F')+Q;
uk1k=rou*(ukk-nz-1)+nz+1;
Uk1k=rou*Ukk;
% Measurement update
xkk=xk1k;
Pkk=Pk1k;
for i=1:N
    % Update IW distribution parameters for Rk
    Bk=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';
    ukk=uk1k+1;
    Ukk=Uk1k+Bk;
    % Update state estimate
    D_R=(ukk-nz-1)\Ukk;
    Pzzk1k=H*Pk1k*H'+D_R;
    Pxzk1k=Pk1k*H';
    Kk=Pxzk1k/Pzzk1k;
    zk1k=H*xk1k;
    xkk=xk1k+Kk*(z-zk1k);
    Pkk=Pk1k-Kk*H*Pk1k;
end
end




