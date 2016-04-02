classdef ReinforcementLearning < handle
    properties
        Action;         %actions or command, for example -5:5
        Action_num;     %total number of actions can be taken
        Action_hist;    %action history (Luca)
        pol_hist;       %policy probabilities history (Luca)

        alpha_v;        %learning rate to update the value function
        alpha_p;        %learning rate to update the policy function
        alpha_n;        %learning rate to update the nature gradient w

        gamma;          %learning rate to update cumulative value or decay rate of moving average
        lambda;         %the regularizatoin factor
        xi;             %discount factor

        Temperature;    %temperature in softmax function in policy network
        weight_range;   %maximum initial weight

        inputDim;       %number of neurons in the input layer
        S1_pol;         %number of neurons in the middle layer of policy network
        S1_val;         %number of neurons in the middle layer of value network
        S2;             %number of actions

        Weights;

        J = 0;          %average of estimated reward
        g;              %"w", the gradient of the policy
        X;              %last input
        Y_pol;          %values of the middle layer in policy network for last input
        Y_val;          %values of the middle layer in value network for last input
        val;            %value estimated for last input
        pol;            %policy for last input
        label_act;      %which action has been chosen
        td;             %TD error
        Weights_hist;   %weights history
        continuous;     %flag whether policy is discrete or continous
    end

    methods
        %PARAM = {Action, alpha_v, alpha_n, alpha_p, xi, gamma, Temperature, lambda, inputDim,
        %         weight_range, loadweights, weights, weights_hist, continuous};
        function obj = ReinforcementLearning(PARAM)
            obj.Action = PARAM{1};

            obj.alpha_v = PARAM{2};
            obj.alpha_n = PARAM{3};
            obj.alpha_p = PARAM{4};
            obj.xi = PARAM{5};
            obj.gamma = PARAM{6};
            obj.Temperature = PARAM{7};
            obj.lambda = PARAM{8};
            obj.inputDim = PARAM{9};
            obj.S2 = length(obj.Action);
            obj.Action_num = length(obj.Action);

            obj.Action_hist = [];
            obj.pol_hist = [];
            obj.td = [];
            obj.Weights_hist = cell(2, 1);
            obj.weight_range = PARAM{10};

            % discrete/continuous action space
            obj.continuous = PARAM{14};

            % load/init
            if (PARAM{11})
                obj.Weights = PARAM{12}(1:2);
                obj.J = PARAM{12}{3};
                obj.g = PARAM{12}{4};
                obj.Weights_hist = PARAM{13};
            else
                obj.NAC_initNetwork();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%initialize the parameters of the class
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function NAC_initNetwork(this)
            % Actor/Policy
            this.Weights{1, 1} = (2 * rand(this.S2, this.inputDim) - 1) * this.weight_range(1);
            % Critic/Value
            this.Weights{2, 1} = (2 * rand(1, this.inputDim) - 1) * this.weight_range(2);

            this.J = 0;                             %average reward (Eq. 3.6, 3.7)
            this.g = zeros(this.S2 * this.inputDim, 1);   %gradient
            this.Weights_hist = this.Weights;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% update the parameters in the network
        %%%
        %%% Xin is the input to the network
        %%% reward is the reward for reinforement learning
        %%% flag_update indicates whether the network should be updated
        %%%
        %%% D contrains the intermedia values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function D = NAC_updateNetwork(this, Xin, reward, flag_update)
            val_new = this.Weights{2, 1} * Xin;
            D = zeros(1, 32);

            if (flag_update) % see Eq. 3.4 to 3.12 in thesis
                this.J = (1 - this.gamma) * this.J + this.gamma * reward;
                delta = reward - this.J + this.xi * val_new - this.val; %TD error
                this.td = delta; %save TD error (Luca)

                dv_val = delta * this.X';
                this.Weights{2, 1} = this.Weights{2, 1} + this.alpha_v * dv_val; %value net update

                dlogv_pol = 1 / this.Temperature * (this.label_act - this.pol) * this.X';
                psi = dlogv_pol(:);

                %alphaforg = this.alpha; %learning rate for w (g)

                deltag = delta * psi - psi * (psi' * this.g);
                this.g = this.g + this.alpha_n * deltag;

                dlambda = this.g;

                dv_pol = reshape(dlambda(1 : numel(dlogv_pol)), size(this.Weights{1, 1}));
                this.Weights{1, 1} = this.Weights{1, 1} * (1 - this.alpha_p * this.lambda);
                this.Weights{1, 1} = this.Weights{1, 1} + this.alpha_p * dv_pol;
                %th2 = 100;
                %L2norm = norm(this.Weights{1, 1},'fro');
                %if (L2norm > th2)
                %    this.Weights{1, 1} = this.Weights{1, 1} * th2/L2norm;
                %end
            end
            this.X = Xin;
            this.val = val_new;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% generate command according to the softmax distribution of the
        %%% output in the policy network
        %%% Xin: feature input to the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function command = softmaxAct(this, Xin)
            poltmp = this.Weights{1, 1} * Xin / this.Temperature;
            this.pol = softmax(poltmp - max(poltmp));   %why?

            assosiatedsoftmax = this.pol;
            softmaxpol = tril(ones(this.Action_num)) * assosiatedsoftmax;
            softmaxpol = softmaxpol / softmaxpol(end) - rand;

            softmaxpol(softmaxpol < 0) = 2;
            [~, index] = min(softmaxpol);
            command = this.Action(index);
            this.label_act = zeros(this.Action_num, 1);
            this.label_act(index) = 1;

            %save action taken and policy values history (Luca)
            this.Action_hist = [this.Action_hist command];
            this.pol_hist = [this.pol_hist this.pol];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pick the action with the maximum possibility amount all the
        %%% commands calculated in policy network
        %%% Xin: feature input to the network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [command, pol, Xin] = Act(this, Xin)
            poltmp = this.Weights{1, 1} * Xin / this.Temperature;
            pol = softmax(poltmp - max(poltmp));
            [~, index] = max(poltmp);
            command = this.Action(index);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Train the reinforcement network for one step
        %%%
        %%% Xin: feature input to the network
        %%% reward is the reward for reinforement learning
        %%% flag_update indicates whether the network should updated
        %%%
        %%% command is the output command
        %%% parameters is the intermedia values keeped for debug
        %%% En is the entropy of policy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [command, paramsC, paramsA] = stepTrain(this, Xin, reward, flag_update)
            this.NAC_updateNetwork(Xin, reward, flag_update);
            command = this.softmaxAct(Xin);
            paramsC = uint8(0);
            paramsA = uint8(0);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the weights during training
        %%% 3rd dim corresponds to iteration, col to weight
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveWeights(this)
            this.Weights_hist{1} = cat(3, this.Weights_hist{1}, this.Weights{1}); %policy net
            this.Weights_hist{2} = cat(3, this.Weights_hist{2}, this.Weights{2}); %value net
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the parameters in a file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveClass(this, configfile)
            weights = cell(2, 1);
            weights{1} = this.Weights;
            weights{2} = this.g;
            save(configfile, 'weights', '-append');
        end
    end
end
