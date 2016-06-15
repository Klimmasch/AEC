%%%
% Natural Actor Critic with Advantage Parameters (discrete RL model)
%%%
classdef ReinforcementLearning < handle
    properties
        actionSpace;    % actions or command, for example -5:5
        nAction;        % total number of admissible actions
        pol_hist;       % policy probabilities history (Luca)

        alpha_v;        % learning rate to update the value function
        alpha_p;        % learning rate to update the policy function
        alpha_n;        % learning rate to update the nature gradient w

        gamma;          % learning rate to update cumulative value or decay rate of moving average
        lambda;         % the regularizatoin factor
        xi;             % discount factor

        temperature;    % temperature in softmax function in policy network
        weight_range;   % maximum initial weight

        input_dim;      % number of neurons in the input layer
        output_dim;     % number of actions

        weightArray;

        J;              % average of estimated reward (Eq. 3.6, 3.7)
        g;              % "w", the gradient of the policy
        feature_prev;   % last input
        val;            % value estimated for last input
        pol;            % policy for last input
        label_act;      % which action has been chosen
        td;             % TD error
        continuous;     % flag whether policy is discrete or continous
    end

    methods
        function obj = ReinforcementLearning(PARAM)
            obj.actionSpace = PARAM{1};

            obj.alpha_v = PARAM{2};
            obj.alpha_n = PARAM{3};
            obj.alpha_p = PARAM{4};
            obj.xi = PARAM{5};
            obj.gamma = PARAM{6};
            obj.temperature = PARAM{7};
            obj.lambda = PARAM{8};
            obj.input_dim = PARAM{9}(1);
            obj.output_dim = length(obj.actionSpace);
            obj.nAction = length(obj.actionSpace);
            obj.weight_range = PARAM{10};
            obj.continuous = PARAM{11};

            % Actor/Policy
            obj.weightArray{1, 1} = (2 * rand(obj.output_dim, obj.input_dim) - 1) * obj.weight_range(1);

            % Critic/Value
            obj.weightArray{2, 1} = (2 * rand(1, obj.input_dim) - 1) * obj.weight_range(1);

            obj.J = 0;
            obj.g = zeros(obj.output_dim * obj.input_dim, 1);
            obj.td = 0;

            obj.pol_hist = [];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% update the parameters in the network
        %%%
        %%% feature is the input to the network
        %%% reward is the reward for reinforement learning
        %%% flag_update indicates whether the network should be updated
        %%%
        %%% D contrains the intermedia values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update(this, feature, reward, flag_update)
            val_new = this.weightArray{2, 1} * feature;

            if (flag_update) % see Eq. 3.4 to 3.12 in thesis
                this.J = (1 - this.gamma) * this.J + this.gamma * reward;
                delta = reward - this.J + this.xi * val_new - this.val; % TD error
                this.td = delta; % save TD error (Luca)

                dv_val = delta * this.feature_prev';
                this.weightArray{2, 1} = this.weightArray{2, 1} + this.alpha_v * dv_val; % value net update

                dlogv_pol = 1 / this.temperature * (this.label_act - this.pol) * this.feature_prev';
                psi = dlogv_pol(:);

                % alphaforg = this.alpha; %learning rate for w (g)

                deltag = delta * psi - psi * (psi' * this.g);
                this.g = this.g + this.alpha_n * deltag;

                dlambda = this.g;

                dv_pol = reshape(dlambda(1 : numel(dlogv_pol)), size(this.weightArray{1, 1}));
                this.weightArray{1, 1} = this.weightArray{1, 1} * (1 - this.alpha_p * this.lambda);
                this.weightArray{1, 1} = this.weightArray{1, 1} + this.alpha_p * dv_pol;
                % th2 = 100;
                % L2norm = norm(this.weightArray{1, 1},'fro');
                % if (L2norm > th2)
                %    this.weightArray{1, 1} = this.weightArray{1, 1} * th2/L2norm;
                % end
            end
            this.feature_prev = feature;
            this.val = val_new;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% generate command according to the softmax distribution of the
        %%% output in the policy network
        %%% feature: feature input to the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function command = softmaxAct(this, feature)
            tmpPol = this.weightArray{1, 1} * feature / this.temperature;
            this.pol = softmax(tmpPol - max(tmpPol)); % why?

            assosiatedSoftmax = this.pol;
            softmaxPol = tril(ones(this.nAction)) * assosiatedSoftmax;
            softmaxPol = softmaxPol / softmaxPol(end) - rand;

            softmaxPol(softmaxPol < 0) = 2;
            [~, index] = min(softmaxPol);
            command = this.actionSpace(index);
            this.label_act = zeros(this.nAction, 1);
            this.label_act(index) = 1;

            % save policy values history (Luca)
            this.pol_hist = [this.pol_hist this.pol];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pick the action with the maximum possibility amount all the
        %%% commands calculated in policy network
        %%% feature: feature input to the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [command, pol, feature] = act(this, feature)
            % tmpPol = this.weightArray{1, 1} * feature / this.temperature;
            tmpPol = this.weightArray{1, 1} * feature / 0.1; % choose a low temperature during testing of policy
            pol = softmax(tmpPol - max(tmpPol));
            [~, index] = max(tmpPol);
            command = this.actionSpace(index);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Train the reinforcement network for one step
        %%%
        %%% feature: feature input to the network
        %%% reward is the reward for reinforement learning
        %%% flag_update indicates whether the network should updated
        %%%
        %%% command is the output command
        %%% parameters is the intermedia values keeped for debug
        %%% En is the entropy of policy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function command = stepTrain(this, feature, reward, flag_update)
            this.update(feature, reward, flag_update);
            command = this.softmaxAct(feature);
        end
    end
end
