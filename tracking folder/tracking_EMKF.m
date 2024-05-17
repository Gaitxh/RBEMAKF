function res = tracking_KF(sim)

H_k = sim.H;
x_ini = sim.inix;
P_ini = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
[x_k_1k_1]=EMKF(F_k,H_k,sim.z,Q_k,R_k,sim.length,x_ini,P_ini);
rec_state = x_k_1k_1;
rec_data = H_k*x_k_1k_1;
res.x = rec_state;
res.z = rec_data;
end

function [zakk]=EMKF(Aa,Ga,data,Qa,Ra,n,inix,iniP)
Cntr = 1;
Lold = 1;
while 1
    Cntr = Cntr + 1;
    %%E-step
    for k = 1:n%forward loop
        yak = data(:,k);
        if k == 1
            zakk1(:,k) = Aa*inix;
            Pzkk1(:,:,k) = Aa*iniP*Aa'+Qa;
        else
            zakk1(:,k) = Aa*zakk(:,k-1);
            Pzkk1(:,:,k) = Aa*Pzkk(:,:,k-1)*Aa'+Qa;
        end
        Kfk(:,:,k) = Pzkk1(:,:,k)*Ga'/(Ra+Ga*Pzkk1(:,:,k)*Ga');
        zakk(:,k) = zakk1(:,k) + Kfk(:,:,k)*(yak-Ga*zakk1(:,k));
        Pzkk(:,:,k) = Pzkk1(:,:,k) - Kfk(:,:,k)*Ga*Pzkk1(:,:,k);
    end
    zakn(:,n) = zakk(:,end);
    Pzkn(:,:,n) = Pzkk(:,:,end);
    for k = n-1:-1:1%backward loop
        Ksk(:,:,k) = Pzkk(:,:,k)*Aa'/Pzkk1(:,:,k+1);
        zakn(:,k) = zakk(:,k)+Ksk(:,:,k)*(zakn(:,k+1)-zakk1(:,k+1));
        Pzkn(:,:,k) = Pzkk(:,:,k) + Ksk(:,:,k)*(Pzkn(:,:,k+1)-Pzkk1(:,:,k+1))*Ksk(:,:,k)';
    end
    Ra = 0*eye(size(yak,1));
    Qa = 0*eye(size(zakk,1));
    Ksk(:,:,n) = Ksk(:,:,end);%%Insufficient data length for data supplementation
    for k = 1:n
        yak = data(:,k);
        Ra = Ra + (yak-Ga*zakn(:,k))*(yak-Ga*zakn(:,k))'+Ga*Pzkn(:,:,k)*Ga';
        if k == 1
            Qa = Qa + (zakn(:,k)-Aa*inix)*(zakn(:,k)-Aa*inix)'...
                +Pzkn(:,:,k)+Aa*iniP*Aa'-Aa*Ksk(:,:,k)*iniP-iniP*Ksk(:,:,k)'*Aa';
        else
            Qa = Qa + (zakn(:,k)-Aa*zakn(:,k-1))*(zakn(:,k)-Aa*zakn(:,k-1))'...
                +Pzkn(:,:,k)+Aa*Pzkn(:,:,k-1)*Aa'-Aa*Ksk(:,:,k)*Pzkn(:,:,k-1)-Pzkn(:,:,k-1)*Ksk(:,:,k)'*Aa';
        end
    end
    L = real(-(n/2)*(sum(log(det(Ra)))+sum(log(det(Qa)))))-.5*(trace(inv(Ra/n)*Ra)+trace(inv(Qa/n)*Qa));
    Ra = Ra/n;
    Qa = Qa/n;
    if Cntr > 20 || (L-Lold)/(L)<1e-4
        break;
    end
    Lold = L;
    Ls(Cntr) = L;
end
end