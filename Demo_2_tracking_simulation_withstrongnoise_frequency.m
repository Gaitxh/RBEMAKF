clc;
clear all;
% close all;
warning off
addpath(genpath(pwd))
% Method_flag = setdiff([1:15],[9,11,12,14]);
Method_flag = [1:9];
Method_Label = {'trueKF','RSTKF','GSTMKF','NVGMNGHVGMKF','SSMKF','AORKF','MLKF', 'LRBEMAKF','STRBEMAKF'};
for title_i = 1:length(Method_flag)
    title_label{title_i} = Method_Label{Method_flag(title_i)};
end
title_label{end+1} = 'Predefine';
Iteration = 1000;
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
intensity = [100,200,300,400,500];
frequency = [.01,.02,.03,.04,.05];
for i = 1:length(frequency)
    sim.fre = frequency(i);
    sim.intensity = intensity(1);
    for Cycle_ii = 1 : Iteration
        sim.inix = [100 100 10 10]';%% defined the initial position and vec
        sim.iniP = diag([100 100 100 100]);%% defined the initial covariance matrix for x0|0
        sim = tracking_simulation_strongnoise(sim);
        sim.inix = sim.inix + utchol(sim.iniP)*randn(sim.m,1);
        for Method_ii = 1:1:length(Method_flag)
            tic
            res = tracking_estimate_outliers(sim,Method_flag(Method_ii));
            [pos{Method_flag(Method_ii)}(Cycle_ii,:),vel{Method_flag(Method_ii)}(Cycle_ii,:)] =...
                indicator_calculate(sim,res);
            T = toc;
            time_spend(Method_flag(Method_ii),Cycle_ii)=T;
            fprintf(['%d,Cycle Time %d, Method %s, Spend Time: %.2f s\n'],...
                i,Cycle_ii,Method_Label{Method_flag(Method_ii)},T);
            %         figure
            %         plot(sim.z(1,:),sim.z(2,:),'m','linewidth',1.5);
            %         hold on
            %         plot(sim.x(1,:),sim.x(2,:),'r','linewidth',2);
            %         hold on
            %         plot(res.x(1,:),res.x(2,:),'black--','linewidth',2);
            %         legend('noised','predefined','estimated');
            %         title(title_label{Method_ii})
        end
    end
    for Method_ii = 1:1:length(Method_flag)
        mean_pos{i}(Method_flag(Method_ii),:)=sqrt(mean(pos{Method_flag(Method_ii)}));
        mean_vel{i}(Method_flag(Method_ii),:)=sqrt(mean(vel{Method_flag(Method_ii)}));
    end
    final_mean_pos(:,i) = mean(mean_pos{i},2);%ARMSE pos
    final_mean_vel(:,i) = mean(mean_vel{i},2);%ARMSE vel
    final_mean_time(:,i) = mean(time_spend,2);%average spend time (second)
end