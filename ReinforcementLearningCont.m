%%%
% Wrapper class for RL actor and critic instances of continuous action space models
%%%
classdef ReinforcementLearningCont < handle
    properties
        weight_range;           % network weight ranges [critic_ji, actor_ji, actor_kj]
        continuous;             % flag whether policy is discrete or continous
        rlFlavour;              % which critic and actor implementation is chosen
        CCritic;                % RL critic instance
        CActor;                 % RL actor instance
        criticLearningRange;    % boundaries for learning rates of critic
        actorLearningRange;     % boundaries for learning rates of actor
        actorDecFac;            % learning rate decay factor for the actor
        criticDecFac;           % learning rate decay factor for the critic
        bias;                   % constant bias in the network
    end

    methods
        function obj = ReinforcementLearningCont(PARAM)
            obj.weight_range = PARAM{10};
            obj.continuous = PARAM{11};
            obj.rlFlavour = PARAM{15};
            obj.criticLearningRange = PARAM{2};
            obj.actorLearningRange = PARAM{4};
            obj.criticDecFac = PARAM{18};
            obj.actorDecFac = PARAM{19};
            obj.bias = PARAM{20};

            % instantiate chosen Actor and Critic
            switch obj.rlFlavour(1)
                case 0
                    %% CACLA
                    % criticParams = {obj.inputDim, obj.weight_range(1), obj.alpha_v, obj.gamma};
                    criticParams = {PARAM{9}(1), obj.weight_range(1), obj.criticLearningRange(1), PARAM{6}};
                    obj.CCritic = CACLACritic(criticParams);
                case 1
                    %% CRG
                    % criticParams = {obj.inputDim, obj.weight_range(1), obj.alpha_v, obj.xi, obj.gamma};
                    criticParams = {PARAM{9}(1), obj.weight_range(1), obj.criticLearningRange(1), PARAM{5}, PARAM{6}};
                    obj.CCritic = CRGCritic(criticParams);
                otherwise
                    error('Critic algorithm [No. #%d] not supported (anymore)!', obj.rlFlavour(1));
            end

            switch obj.rlFlavour(2)
                case 0
                    %% CACLAVar [Lukas's interpretation of CACLA appoach]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec, obj.regularizer};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{7}, PARAM{12}, PARAM{13}, PARAM{16}, PARAM{17}};
                    obj.CActor = CACLAVarActorLu(actorParams);
                case 1
                    %% CACLAVar [Alex's interpretation of CACLA appoach]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{7}, PARAM{12}, PARAM{13}, PARAM{16}};
                    obj.CActor = CACLAVarActorAl(actorParams);
                case 2
                    %% CACLAVar [CACLA appoach with std. Backpropagation]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{7}, PARAM{12}, PARAM{13}, PARAM{16}};
                    obj.CActor = CACLAVarActorBp(actorParams);
                case 3
                    %% CACLAVar
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{7}, PARAM{12}, PARAM{13}, PARAM{16}};
                    obj.CActor = CACLAVarActor(actorParams);
                case 4
                    %% CACLAVar2
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{7}, PARAM{12}, PARAM{13}, PARAM{16}};
                    obj.CActor = CACLAVarActor2(actorParams);
                % case 5
                %     % TODO: unsupported yet
                %     %% CNGFI
                %     % actorParams = {obj.inputDim, obj.outputDim, obj.weight_range(2:3), obj.alpha_p, obj.alpha_v, obj.varianceRange, obj.fiScale, obj.varDec};
                %     actorParams = {PARAM{9}, 1, obj.weight_range(2:3), obj.actorLearningRange(1), PARAM{2}, PARAM{7}, PARAM{14}, PARAM{16}};
                %     obj.CActor = CNGFIActor(actorParams);
                otherwise
                    error('Actor algorithm [No. #%d] not supported (anymore)!', obj.rlFlavour(2));
            end
        end

        %%% Generate a pure action output without Gaussean policy noise
        %   @param feature: input feature vector
        function command = act(this, feature)
            command = this.CActor.actHard(feature);
        end

        %%% Train the critic's and actor's networks and generate action output according to policy
        %   @param feature:         input feature vector
        %   @param reward:          reward signal from sparse coders
        %   @param flag_update:     indicates whether weights shall be updated (additional fine grained control)
        %   @return:                delta muscle, i.e. change in eye muscle excitation(s)
        function command = stepTrain(this, feature, reward, flag_update)
            this.CCritic.train(feature, reward, flag_update);
            command = this.CActor.train(feature, this.CCritic.delta, flag_update);
        end
    end
end
