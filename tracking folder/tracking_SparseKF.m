function res = tracking_KF(sim)

H_k = sim.H;
x_k_1k_1 = sim.inix;
P_k_1k_1 = sim.iniP;
F_k = sim.F;
Q_k = sim.Q0;
R_k = sim.R0;
for k = 1 : sim.length
    z_k = sim.z(:,k);
    [x_k_1k_1,P_k_1k_1]=SKF(z_k,H_k,F_k,x_k_1k_1,P_k_1k_1,Q_k,R_k,size(F_k,1),size(H_k,1));
    rec_state(:,k) = x_k_1k_1;
    rec_data(:,k) = H_k*x_k_1k_1;
end
res.x = rec_state;
res.z = rec_data;
end

function [x_k,P_kx] = ...
    SKF(y_k,C,A,x_k_1k_1,Pk_1x,Q,R,m,n)
%% 参数预定义
B = eye(m);
A1= eye(m);

regre = C*B;
x_kk_1 = A*x_k_1k_1;%---(9)
Sigma = (C*A)*Pk_1x*(C*A)'+R;%---(31)
M_k = (A1+regre'/Sigma*regre)\regre'/Sigma;%---(38)
f_k_1 = M_k*(y_k-C*x_kk_1);%---(10)
x_kkx = x_kk_1+B*f_k_1;%---(11)

%% 更新相应的alpha值，实现状态的稀疏
t = y_k-C*x_kk_1;%---(30)
regre = C*B;%---(30)
[U,H,V] = svd(Sigma);%SVD分解
T = inv(H.^(1/2))*U;%---(33)
regre_1 = T*regre;%---(33)
t_1 = T*t;%---(33)
alpha = alpha_calculate(regre_1,t_1,eye(n),m);%algorithm-1 计算alpha的超参值
A1 = diag(alpha);

%% 更新相应参数
miu = (A1+regre'/Sigma*regre)\regre'/Sigma*t;%---(37)
M_k = (A1+regre'/Sigma*regre)\regre'/Sigma;%---(38)
I = eye(m);
P_k_1f = inv(A1);%---(42)
P_kxx = (I-B*M_k*C)*A*Pk_1x*A'*(I-B*M_k*C)'...
    +B*(I-M_k*regre)*P_k_1f*(I-M_k*regre)'*B'+B*M_k*R*(B*M_k)';%---(41)
S_kx = -B*M_k*R;%---(47)
R_kx = C*P_kxx*C'+R+C*S_kx+S_kx'*C';%---(45)
V_kx = P_kxx*C'+S_kx;%---(46)
K_k = V_kx/R_kx;%---(48)
P_kx = K_k*R_kx*K_k'-V_kx*K_k'-K_k*V_kx'+P_kxx;%---(44)
x_k = x_kkx+K_k*(y_k-C*x_kkx);%---(12)
end

function alpha = alpha_calculate(regre,t,Sigma,m)
for i = 1:m
    res_i = setdiff(1:m,i);
    alpha = inf*ones(m,1);%所有其他的alpha都被设置为无穷大
    alpha(i) = norm(regre(:,i),2)/(norm(regre(:,i)'*t,2)/norm(regre(:,i),2)-1);%计算初值
    C_i = (Sigma+regre(:,res_i)*diag(1./alpha(res_i))*regre(:,res_i)');
    s(i) = regre(:,i)'/C_i*regre(:,i);
    q(i) = regre(:,i)'/C_i*t;
    theta = q(i).^2-s(i);
    if theta> 0
        alpha_ensure(1,i) = s(i).^2/theta;
    elseif theta <= 0
        alpha_ensure(1,i) = inf;
    end
end
count = 1;
while 1
    count = count + 1;
    for i = 1:m
        res_i = setdiff(1:m,i);
        C_i = (Sigma+regre(:,res_i)*diag(1./alpha_ensure(count-1,res_i))*regre(:,res_i)');
        s(i) = regre(:,i)'/C_i*regre(:,i);
        q(i) = regre(:,i)'/C_i*t;
        theta = q(i).^2-s(i);
        if theta > 0 && alpha_ensure(count-1,i) < inf
            alpha_ensure(count,i) = s(i).^2/theta;
        elseif theta > 0 && alpha_ensure(count-1,i) == inf
            alpha_ensure(count,i) = inf;
        elseif theta <= 0
            alpha_ensure(count,i) = inf;
        end
    end
    if count > 10
        break;
    end
end
alpha = alpha_ensure(count,:);
% if length(alpha(find(alpha~=inf)))~=0
%     alpha
% end
end