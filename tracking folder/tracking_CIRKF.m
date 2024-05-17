function res = tracking_CIRKF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    f_v=5.0;N = 10;
    [x_k_1k_1,P_k_1k_1]=CIRKF(x_k_1k_1,P_k_1k_1,z_k,H_k,F_k,Q_k,R_k,f_v,N);
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end
function [xkk,Skk]=CIRKF(xkk,Skk,z,H,F,Q,R,v,N)
xold = xkk;
%%%%%%Time update
nx=size(xkk,1);    

nz=size(z,1);

nPts=2*nx;         

Xk1k1=CR(xkk,Skk);

Xkk1=ckf_ProssEq(Xk1k1,F);                                

xkk1=sum(Xkk1,2)/nPts;       

Pkk1=Xkk1*Xkk1'/nPts-xkk1*xkk1'+Q;             

%%%%%%Measurement update
lamda=1;

for i=1:N
    
    %%%%%%Update state
    Xi=CR(xkk1,Pkk1);
    
    Zi=ckf_Mst(Xi,nz,H);
    
    zkk1=sum(Zi,2)/nPts;  
    
    Pzzkk1=Zi*Zi'/nPts-zkk1*zkk1'+R/lamda;   
    
    Pxzkk1=Xi*Zi'/nPts-xkk1*zkk1';

    Wk=Pxzkk1*inv(Pzzkk1);      

    xkk=xkk1+Wk*(z-zkk1);
    
    Skk=Pkk1-Wk*Pzzkk1*Wk';   
    
    %%%%%%Update lamda
    Xkk=CR(xkk,Skk);
    
    Zkk=ckf_Mst(Xkk,nz,H);
    
    Zkk=repmat(z,1,nPts)-Zkk;
    
    Dk=Zkk*Zkk'/nPts;
    
    gama=trace(Dk*inv(R));
         
    lamda=(v+nz)/(v+gama);
    
end

end
