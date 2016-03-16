classdef ReinforcementLearningCont < handle
    properties
        alpha_v;        %learning rate to update the value function
        alpha_p;        %learning rate to update the policy function
        alpha_n;        %learning rate to update the nature gradient w

        gamma;          %learning rate to update cumulative value or decay rate of moving average
        lambda;         %the regularizatoin factor
        xi;             %discount factor

        variance;    %temperature in softmax function in policy network
        weight_range;   %maximum initial weight
        S0;             %number of neurons in the input layer

        % DEPRECATED
        % TODO: update model reload
        % Weights;
        % J = 0;          %average of estimated reward
        % g;              %intermedia variable to keep track of "w" which is the gradient of the policy
        % Weights_hist;   %weights history

        continuous;     %flag whether policy is discrete or continous

        % for continuous action space with gaussian policy
        CCritic;
        CActor;
    end

    methods
        function obj = ReinforcementLearningCont(PARAM)
            obj.alpha_v = PARAM{2};
            obj.alpha_n = PARAM{3};
            obj.alpha_p = PARAM{4};
            obj.xi = PARAM{5}; % in continuous case used as gamma
            obj.gamma = PARAM{6}; % in continuous case used as eta
            obj.variance = PARAM{7};
            obj.lambda = PARAM{8};
            obj.S0 = PARAM{9};
            obj.weight_range = PARAM{10};

            % continuous action space (lukas)
            obj.continuous = PARAM{14};

            %CriticParams = {alpha_v, eta, gamma, featureDimension, initialWeightRange}
            %CriticParams = {0.01, 0.05, 0.3, 576, 0.15}; --> params from original implementation
            % CriticParams = {obj.alpha_v, obj.gamma, obj.xi, obj.S0, obj.weight_range(1)};
            % obj.CCritic = CCriticG(CriticParams);
            CriticParams = {obj.S0, obj.weight_range(1), obj.alpha_v, obj.xi, obj.gamma};
            obj.CCritic = CRGCritic(CriticParams);

            %ActorParams = {alpha_p, alpha_n, featureDimension, initialWeightRange, actorHiddenType, variance};
            %ActorParams = {0.001, 0.01, 576, 0.22, 'tanh'}; original params
            % ActorParams = {obj.alpha_p, obj.alpha_n, obj.S0, obj.weight_range(2:3), 'tanh', 'default', obj.variance};
            % obj.CActor = CActorG(ActorParams);
            ActorParams = {obj.S0, 1, obj.weight_range(2:3), obj.alpha_p, obj.variance};
            obj.CActor = CRGActor(ActorParams);

            % DEPRECATED
            % TODO: update model reload
            % obj.Weights_hist = cell(2, 1);
            % load/init
            if (PARAM{11})
                sprintf('Model reloading function is DEPRECATED and therefore currently not supported!');
                return;
                % obj.Weights = PARAM{12}(1:2);
                % obj.J = PARAM{12}{3};
                % obj.g = PARAM{12}{4};
                % obj.Weights_hist = PARAM{13};
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% generate command according to the softmax distribution of the
        %%% output in the policy network
        %%% Xin: feature input to the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function command = softmaxAct(this, Xin)
            command = this.CActor.actHard(Xin);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Train the reinforcement network for one step
        %%%
        %%% feature is the input to the network
        %%% reward is the reward for reinforement learning
        %%% flag_update indicates whether the network should updated
        %%%
        %%% command is the output command
        %%% parameters is the intermedia values keeped for debug
        %%% En is the entropy of policy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function command = stepTrain(this, feature, reward, flag_update)
            this.CCritic.train(feature, reward, flag_update);
            command = this.CActor.train(feature, this.CCritic.delta, flag_update);
        end

        % DEPRECATED
        % TODO: update model reload
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the weights during training
        %%% 3rd dim corresponds to iteration, col to weight
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveWeights(this)
            sprintf('Model reloading function is DEPRECATED and therefore currently not supported!');
            return;
            % this.Weights_hist{1} = cat(3, this.Weights_hist{1}, this.Weights{1}); %policy net
            % this.Weights_hist{2} = cat(3, this.Weights_hist{2}, this.Weights{2}); %value net
        end

        % DEPRECATED
        % TODO: update model reload
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the parameters in a file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveClass(this, configfile)
            sprintf('Model reloading function is DEPRECATED and therefore currently not supported!');
            return;
            % weights = cell(2, 1);
            % weights{1} = this.Weights;
            % weights{2} = this.g;
            % save(configfile, 'weights', '-append');
        end
    end
end
