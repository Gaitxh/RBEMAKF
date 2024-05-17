function sim = tracking_simulation_weeknoise(sim)
xk_1 = sim.inix;
for t = 1:sim.length
    %%%%True noise covariance matrices
    Q=(6.5+0.5*cos(pi*t/sim.length))*sim.Q1;
    R=(0.1+0.05*cos(pi*t/sim.length))*sim.R1;
    sim.Q(:,:,t) = Q;
    sim.R(:,:,t) = R;
    %%%%Square-root of noise covariance matrices
    SQ=utchol(Q);
    SR=utchol(R);
    xk = sim.F*xk_1+SQ*randn(sim.m,1);
    zk = sim.H*xk+SR*randn(sim.n,1);
    xk_1 = xk;
    
    %% save the initial state value
    sim.x(:,t) = xk;
    sim.z(:,t) = zk;
end
% figure
% plot(sim.z(1,:),sim.z(2,:))
% hold on
% plot(sim.x(1,:),sim.x(2,:))
% legend('observation','ture')
end