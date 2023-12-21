clc;
clear all;
% close all;
warning off
addpath(genpath(pwd))
% Method_flag = setdiff([1:15],[9,11,12,14]);
Method_flag = [1:8];
Method_Label = {'trueKF','RCKF','VBAKF','UPAKF','SBAKF','EMAKF','SKF','RBEMAKF'};
for title_i = 1:length(Method_flag)
    title_label{title_i} = Method_Label{Method_flag(title_i)};
end
title_label{end+1} = 'Predefine';
Iteration = 100;
sim.length = 200;%% the length of time series
%% the parameters that prefined for simulation data generation
sim.n = 2;sim.m = 4;%% defined the dismension of the dynamic system
sim.H = [eye(sim.n) zeros(sim.n)];
sim.F = [eye(2) eye(2);zeros(2) eye(2)];%% defined the trasform matrix of the processment formulation
%% defined the initial processment and measurement noise covariance matrices
T=1;q=1;r=100;
sim.Q1=[T^3/3*eye(sim.n) 0*eye(sim.n);0*eye(sim.n) T*eye(sim.n)]*q;sim.R1=r*[1 0;0 1];
%% defined the initial input processment and measurement noises covariance matrices
alfa=1;beta=100;
sim.Q0=alfa*eye(sim.m);sim.R0=beta*eye(sim.n);
%% the iteration of the KF comparison
for Cycle_ii = 1 : Iteration
    sim.inix = [100 100 10 10]';%% defined the initial position and vec
    sim.iniP = diag([100 100 100 100]);%% defined the initial covariance matrix for x0|0
    sim = tracking_simulation_weeknoise(sim);
    sim.inix = sim.inix + utchol(sim.iniP)*randn(sim.m,1);
    for Method_ii = 1:1:length(Method_flag)
        tic
        res = tracking_estimate(sim,Method_flag(Method_ii));
        [pos{Method_flag(Method_ii)}(Cycle_ii,:),vel{Method_flag(Method_ii)}(Cycle_ii,:)] =...
            indicator_calculate(sim,res);
        T = toc;
        time_spend(Method_flag(Method_ii),Cycle_ii)=T;
        fprintf(['Cycle Time %d, Method %s, Spend Time: %.2f s\n'],...
            Cycle_ii,Method_Label{Method_flag(Method_ii)},T);
    end
end
for Method_ii = 1:1:length(Method_flag)
    mean_pos(Method_flag(Method_ii),:)=sqrt(mean(pos{Method_flag(Method_ii)}));
    mean_vel(Method_flag(Method_ii),:)=sqrt(mean(vel{Method_flag(Method_ii)}));
end
final_mean(:,1) = mean(mean_pos,2);%ARMSE pos
final_mean(:,2) = mean(mean_vel,2);%ARMSE vel
final_mean(:,3) = mean(time_spend,2);%average spend time (second)


colors = [
    0.0863, 0.3137, 0.5451; % 深蓝绿色
    0.1529, 0.5412, 0.4431; % 深绿色
    0.3765, 0.6235, 0.4824; % 浅绿色
    0.4980, 0.6863, 0.6157; % 淡蓝绿色
    0.7922, 0.8118, 0.5882; % 浅黄绿色
    0.9294, 0.6941, 0.1255; % 橙色
    0.9608, 0.4863, 0.2275; % 橙红色
    0.7373, 0.1961, 0.1765; % 深红色
    0.6627, 0.3294, 0.3137; % 更深红色 (突出显示)
    0.6157, 0.1725, 0.3216; % 深紫色
    0.4667, 0.2588, 0.4314; % 紫色
    0.3216, 0.3216, 0.3216; % 灰色
];
H1 = figure
subplot(2,2,1)
for Method_ii = 1:1:7
    plot(mean_pos(Method_flag(Method_ii),:),'Color', colors(Method_ii, :),'linewidth',2);
    hold on
end
plot(mean_pos(8,:),'r','linewidth',2);
legend(title_label);
title('RMSE pos')
axis([1 200 2.5 7])
subplot(2,2,3)
for Method_ii = 1:1:7
    plot(mean_vel(Method_flag(Method_ii),:),'Color', colors(Method_ii, :),'linewidth',2);
    hold on
end
plot(mean_vel(8,:),'r','linewidth',2);
title('RMSE vel')
axis([1 200 4 10])
width = 1000;  % 800 像素宽度
height = 1000; % 600 像素高度
set(H1, 'Position', [100, 100, width, height]);

H1 = figure
subplot(6,2,1)
for Method_ii = 1:1:7
    plot(mean_pos(Method_flag(Method_ii),:),'Color', colors(Method_ii, :),'linewidth',2);
    hold on
end
plot(mean_pos(8,:),'r','linewidth',2);
legend(title_label);
title('RMSE pos')
subplot(6,2,3)
for Method_ii = 1:1:7
    plot(mean_vel(Method_flag(Method_ii),:),'Color', colors(Method_ii, :),'linewidth',2);
    hold on
end
plot(mean_vel(8,:),'r','linewidth',2);
title('RMSE vel')
width = 1000;  % 800 像素宽度
height = 1000; % 600 像素高度
set(H1, 'Position', [100, 100, width, height]);
