%%%
% Continuous Actor Critic Learning Automaton Critic
%%%
classdef CACLACritic < handle
    properties
        % network parameters
        input_dim;
        v_ji;           % input -> output weights
        v_init_range;   % weight initialization
        regularizer;    % regularization of weights

        % Reinforcement Learning parameters
        delta;          % estimate of Temporal Difference error
        alpha_v;        % step-size schedule of v_ji update "learning rate"
        gamma;          % discount factor of previous state
        value;          % critic's value estimation (output)

        % input/output tracking of previous time step
        feature_prev;   % feature vector

        % model history tracking
        params;
    end

    methods
        function obj = CACLACritic(PARAM)
            obj.input_dim = PARAM{1};
            obj.v_init_range = PARAM{2};
            obj.alpha_v = PARAM{3};
            obj.gamma = PARAM{4};
            obj.regularizer = PARAM{5};

            obj.v_ji = (2 * rand(1, obj.input_dim) - 1) * obj.v_init_range; % [-1, 1] * v_init_range
            % obj.v_ji = rand(1, obj.input_dim) * obj.v_init_range; % [0, 1] * v_init_range

            obj.delta = 0;
            obj.feature_prev = zeros(obj.input_dim, 1);
            obj.params = zeros(1, 2);
        end

        function calculate(this, feature)
            this.value = this.v_ji * feature;
        end

        function update(this, reward)
            this.delta = reward + this.gamma * this.value - this.v_ji * this.feature_prev;
            % this.delta = (reward * (1 - this.gamma)) + this.gamma * this.value - this.v_ji * this.feature_prev; % scale the reward to make the value indep. of gamma
            dv_ji = this.alpha_v * this.delta * this.feature_prev';
            this.v_ji = (1 - (this.regularizer * this.alpha_v)) * this.v_ji; % similar to the actor 
            this.v_ji = this.v_ji + dv_ji;
        end

        function train(this, feature, reward, flag_update)
            this.calculate(feature);
            if (flag_update)
                this.update(reward);
            end
            this.feature_prev = feature;
            this.params(1) = sum(sum(abs(this.v_ji)));
            this.params(2) = sum(sum(this.v_ji .^ 2));
        end
    end
end
