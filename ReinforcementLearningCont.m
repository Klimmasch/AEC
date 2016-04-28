classdef ReinforcementLearningCont < handle
    properties
        weight_range;
        continuous;     %flag whether policy is discrete or continous
        rlFlavour;      %which critic and actor implementation is chosen
        CCritic;
        CActor;
    end

    methods
        function obj = ReinforcementLearningCont(PARAM)
            obj.weight_range = PARAM{10};
            obj.continuous = PARAM{14};
            obj.rlFlavour = PARAM{18};

            % instantiate chosen Actor and Critic
            switch obj.rlFlavour(1)
                case 0
                    %% CACLA
                    % criticParams = {obj.inputDim, obj.weight_range(1), obj.alpha_v, obj.gamma};
                    criticParams = {PARAM{9}(1), obj.weight_range(1), PARAM{2}, PARAM{6}};
                    obj.CCritic = CACLACritic(criticParams);
                case 1
                    %% CRG
                    % criticParams = {obj.inputDim, obj.weight_range(1), obj.alpha_v, obj.xi, obj.gamma};
                    criticParams = {PARAM{9}(1), obj.weight_range(1), PARAM{2}, PARAM{5}, PARAM{6}};
                    obj.CCritic = CRGCritic(criticParams);
                otherwise
                    sprintf('Critic algorithm not supported (anymore)!')
                    return;
            end

            switch obj.rlFlavour(2)
                case 0
                    %% CACLAVar [Lukas's interpretation of CACLA appoach]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), PARAM{4}, PARAM{7}, PARAM{15}, PARAM{16}, PARAM{19}};
                    obj.CActor = CACLAVarActorLu(actorParams);
                case 1
                    %% CACLAVar [Alex's interpretation of CACLA appoach]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), PARAM{4}, PARAM{7}, PARAM{15}, PARAM{16}, PARAM{19}};
                    obj.CActor = CACLAVarActorAl(actorParams);
                case 2
                    %% CACLAVar [CACLA appoach with std. Backpropagation]
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), PARAM{4}, PARAM{7}, PARAM{15}, PARAM{16}, PARAM{19}};
                    obj.CActor = CACLAVarActorBp(actorParams);
                case 3
                    %% CACLAVar
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), PARAM{4}, PARAM{7}, PARAM{15}, PARAM{16}, PARAM{19}};
                    obj.CActor = CACLAVarActor(actorParams);
                case 4
                    %% CACLAVar2
                    % actorParams = {[obj.inputDim, obj.hiddenDim, obj.outputDim], obj.weight_range(2:3), obj.alpha_p, obj.varianceRange, obj.deltaVar, obj.eta, obj.varDec};
                    actorParams = {PARAM{9}, obj.weight_range(2:3), PARAM{4}, PARAM{7}, PARAM{15}, PARAM{16}, PARAM{19}};
                    obj.CActor = CACLAVarActor2(actorParams);
                % case 5
                %     % TODO: unsupported yet
                %     %% CNGFI
                %     % actorParams = {obj.inputDim, obj.outputDim, obj.weight_range(2:3), obj.alpha_p, obj.alpha_v, obj.varianceRange, obj.fiScale, obj.varDec};
                %     actorParams = {PARAM{9}, 1, obj.weight_range(2:3), PARAM{4}, PARAM{2}, PARAM{7}, PARAM{17}, PARAM{19}};
                %     obj.CActor = CNGFIActor(actorParams);
                otherwise
                    sprintf('Actor algorithm not supported (anymore)!')
                    return;
            end

            % DEPRECATED
            % TODO: update model reload
            % obj.Weights_hist = cell(2, 1);
            % load/init
            if (PARAM{11})
                sprintf('Model reloading function is DEPRECATED and therefore currently not supported!')
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
        function command = act(this, Xin)
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
            sprintf('Model reloading function is DEPRECATED and therefore currently not supported!')
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
            sprintf('Model reloading function is DEPRECATED and therefore currently not supported!')
            return;
            % weights = cell(2, 1);
            % weights{1} = this.Weights;
            % weights{2} = this.g;
            % save(configfile, 'weights', '-append');
        end
    end
end
