function GPmain()

clear all; clc
dimension_Feature = 41;
initialWeightRange = [0.15,0.22];
%  value network
alpha_v = 0.01;
eta = 0.05;
gamma = 0.3; %0.8;
% PARAMCritic = {alpha_v,eta,gamma,dimension_Feature,initialWeightRange(1)};
% critic = CCritic(PARAMCritic);
PARAMCritic = {dimension_Feature, 0.0244, 0.9, 0.3, 0.3};
critic = CRGCritic(PARAMCritic);
% policy network
alpha_n = 0.01;
alpha_p = 0.001;
actor_hidden_type = 'tanh';
% PARAMActor = {alpha_p,alpha_n,dimension_Feature,initialWeightRange(2),actor_hidden_type};
% actor = CActorG(PARAMActor);
PARAMActor = {dimension_Feature, 1, [0.0244, 0.0244], 0.8, 0.5};
actor = CRGActor(PARAMActor);
clearvars -except critic actor;

%% training
rng(23);
T_train = 100000;
Vd_train = zeros(1,T_train);

V0 = 0;
minV = -8;
maxV = 8;
% initialization of trajectory
V_train = zeros(1,T_train);
V_train(1) = V0;
act = 0;

% weights_track = zeros(8,T_train);
weights_track = zeros(2, T_train);
d_track = zeros(1,T_train);
delta_track = zeros(1,T_train);
reward_track = zeros(1,T_train);
criticZ_track = zeros(1,T_train);
features = 0.1*ones(41,41) + eye(41);

verbose = 0;

degrees = load('Degrees.mat');              %loads tabular for resulting degrees as 'results_deg'

command = 0;
command_prev = 0;

%%% Helper function that maps muscle activities to resulting angle
function [angle] = getAngle(command)
    cmd = (command * 10) + 1;                               % scale commands to table entries
    angle = interp2(degrees.results_deg, 1, cmd);   % interpolate in tabular
end

function [cmd] = checkCmd(cmd)
    i0 = cmd < 0;
    cmd(i0) = 0;
    i1 = cmd > 1;
    cmd(i1) = 1;
end

% training
for t=1:T_train
    if(mod(t-1,50)==0)
        x = randi([-5 5],1);
%         x = (2*rand()-1);
    end
    Vd_train(t) = x; %vergence desired
    if (verbose)
        display(x)
    end
    % Vd_train(t) = 15*sin(t/30);
    % generate action
    if(t>1)
        command = getAngle(checkCmd(command_prev + act));
%         command = minV + (maxV - minV) * rand;
        V_train(t) = command;
        % V_train(t) = V_train(t-1) + act;
        if (verbose)
            display(act)
        end
    end

    % if ( (V_train(t)>maxV) || (V_train(t)<minV) )
    if (V_train(t) > maxV)
        % V_train(t) = randi([-5 5],1);
        V_train(t) = maxV;
    elseif (V_train(t) < minV)
        V_train(t) = minV;
    end
    if (verbose)
        sprintf('V_train(%d) = %g', t, V_train(t))
    end

    % generate vector and reward
    d = Vd_train(t) - V_train(t); % actual disparity d
%     d = round(Vd_train(t) - V_train(t)); % actual disparity d
    Xin = features(:,round(d)+21); % feature vector, d+21th column
%     Reward = exp(-d^2/20); % the less the disparity the greater the reward, i.e. the less the punishment
    Reward = exp(-d^2/20);% - 0.001* critic.params(1) - 0.001 * actor.params(1);
    if (verbose)
        display(Reward)
    end
%     [delta,param_c] = critic.train(Xin,Reward,t);
    critic.train(Xin,Reward,mod(t-1,50));
%     [act,param_a] = actor.train(Xin,delta,t);
    act = actor.train(Xin,critic.delta,mod(t-1,50));
%     act = minV + (maxV - minV) * rand;
    %
%     weights_track(:,t) = [param_c,param_a]';
    weights_track(:,t) = [critic.params(1), actor.params(1)];
    d_track(t) = d;
    delta_track(t) = critic.delta;
    reward_track(t) = Reward;
%     criticZ_track(t) = param_c(2);
    if(~mod(t,ceil(T_train/100)))
        disp([num2str(t*100/T_train) '% is finished']);
    end
%      plot(delta_track);
%     plot(reward_track);
% t
end
windowSize = T_train/100;
ind = 1:T_train;
vergerr = filter(ones(1, windowSize) / windowSize, 1, abs(d_track(ind)));
figure; grid on; grid minor; plot(vergerr); title('moving average of disparity'); xlabel('time');ylabel('disparity');
figure; plot([Vd_train(:);V_train(:)]'); title('vergence of agent vs desired vergence');
% figure; plot(d_track(:)); title('disparity');
% figure; plot(delta_track(:)); title('delta');
% figure; plot(reward_track(:)); title('reward');
%correlations of critic: delta positively correlates with the reward, Z is
%learned to be correlated to the reward
%figure; scatter3(delta_track, reward_track, criticZ_track, 9, 1:T_train); xlabel('delta');ylabel('reward');zlabel('Z');
% NAC_NN.NAC_saveWeights();

% figure(1)
% subplot(2,1,1)
% plot([1:T_train], Vd_train);
% axis([1 T_train min(Vd_train) max(Vd_train)]);
% subplot(2,1,2)
% plot([1:T_train], V_train);
% axis([1 T_train min(V_train) max(V_train)]);

% figure(2)
% plot([1:T_train], [AMem;deltaMem]);
% axis([1 T_train min(Vd_train) max(Vd_train)]);

% figure(3)
% subplot(3,1,1)
% plot([1:T_train], valueMem);
% axis([1 T_train min(valueMem) max(valueMem)]);
% subplot(3,1,2)
% plot([1:T_train], JMem);
% axis([1 T_train min(JMem) max(JMem)]);
% subplot(3,1,3)
% plot([1:T_train], deltaMem);
% axis([1 T_train min(deltaMem) max(deltaMem)]);

%% testing

% flag_train = 0;
%
% T_test = 500;
% % initialization of trajectory
% V1_test = zeros(1,T_test);
% V1_test(1) = V0;
%
% % training
% for t = 1:T_test
%     % generate action
%     if (t>1)
%         act = NAC_NN.NAC_genAction(flag_train);
%         V1_test(t) = V1_test(t-1)+act;
%     else
%         V1_test(t) = V0;
%     end
%     % generate input vector and reward
%     % *** APPROXIMATED NOW, NEED TO CHANGE TO MODEL OR SPARSE CODING ***
%     d = Vd_train(t) - V1_test(t);
%     Xin = d;
%     % update NAC network
%     NAC_NN.NAC_updatePolicy(Xin);
% end
%
% % plot
% figure(2);
% plot(1:T_test, [Vd_train(1:T_test);V1_test]);
%plot(1:T_test,V1_test);

% function plotCritic
    disps = minV:1:maxV;
    values = zeros(size(disps,2),1);
    commands = zeros(size(disps,2),1);
    for d = disps
        xin = features(:,d+21);
        critic.calculate(xin);
        values(d+abs(minV)+1) = critic.value;
        commands(d+abs(minV)+1) = actor.act(xin);
    end
    figure;
    plot(disps,values);
    xlabel('disparities');
    ylabel('estimated value');
    grid on;
    figure;
    plot(disps,commands);
    xlabel('disparities');
    ylabel('relative command');
    grid on;
% end

% plotCritic

end


