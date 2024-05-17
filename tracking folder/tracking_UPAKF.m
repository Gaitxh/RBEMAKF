function res = tracking_KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
alfa=5;
beta=5;
Nm=50;
rou=1-exp(-5);
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1,alfa,beta]=akf_upml(x_k_1k_1,P_k_1k_1,F_k,H_k,z_k,Q_k,R_k,alfa,beta,rou,Nm);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [xk1k1,Pk1k1,alfa,beta] = akf_upml(xkk,Pkk,F,H,y,Q,R,alfa0,beta0,rou,N)

%%%%%%%%% Time update
xk1k=F*xkk;
Pk1k=F*Pkk*F'+Q;

%%%%%%%%% Initialization
xk1k1=xk1k;
Pk1k1=Pk1k;

alfa1=rou*alfa0;
beta1=rou*beta0;

E_gamma=beta1/(alfa1+beta1);

E_log_tau=psi(alfa1)-psi(alfa1+beta1);
E_log_1_tau=psi(beta1)-psi(alfa1+beta1);

%%%%%%%%% Variational measurement update
for i=1:N
    
    %%%%%%%%%%
    xk1k1_i=xk1k1;
    
    %%%%%%%%% Update state vector
    yk1k=H*xk1k;
    
    D_R=R/E_gamma;
    
    Sk=H*Pk1k*H';
    
    Pyyk1k=Sk+D_R;
    
    Pxyk1k=Pk1k*H';
    
    Kk=Pxyk1k*inv(Pyyk1k);
    
    xk1k1=xk1k+Kk*(y-yk1k);
    
    Pk1k1=Pk1k-Kk*Pyyk1k*Kk';
    
    %%%%%%%%% Determine convergence
    td=norm(xk1k1-xk1k1_i)/norm(xk1k1);
    
    if td<=1e-16
        break;
    end
    
    %%%%%%%%% Update gamma
    Ak=(y-H*xk1k1)*(y-H*xk1k1)'+H*Pk1k1*H';
    
    Bk=y*y';
    
    Pr1=exp(E_log_1_tau-0.5*trace(Ak*inv(R)));
    
    Pr2=exp(E_log_tau-0.5*trace(Bk*inv(R)));
    
    %%%%%%%%% Update tau
    alfa=alfa1-E_gamma+1;
    
    beta=beta1+E_gamma;
    
    %%%%%%%%% Compute the expectations
    E_gamma=Pr1/(Pr1+Pr2);
    
    E_log_tau=psi(alfa)-psi(alfa+beta);
    
    E_log_1_tau=psi(beta)-psi(alfa+beta);
    
end
end