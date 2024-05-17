function sim = tracking_simulation(sim)
xk_1 = sim.inix;
Q1 = sim.Q1;R1 = sim.R1;
for t = 1:sim.length
    %%%%True noise covariance matrices
    Q=(6.5+0.5*cos(pi*t/sim.length))*Q1;
    R=(0.1+0.05*cos(pi*t/sim.length))*R1;
    sim.Q(:,:,t) = Q;
    sim.R(:,:,t) = R;
    %%%%Square-root of noise covariance matrices
    SQ=utchol(Q);
    SR=utchol(R);
    xk = sim.F*xk_1+utchol(SQ)*randn(sim.m,1);
    zk = sim.H*xk+utchol(SR)*randn(sim.n,1);
    xk_1 = xk;
    
    %% save the initial state value
    sim.x(:,t) = xk;
    sim.z(:,t) = zk;
end
% plot(sim.z(1,:),sim.z(2,:))
% hold on
% plot(sim.x(1,:),sim.x(2,:))
% legend('observation','ture')
end