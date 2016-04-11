%%% Continuous Actor Critic Learning Automaton Actor
classdef CACLAActor < handle
    properties
        % network parameters
        input_dim;
        hidden_dim;
        output_dim;
        wp_ji;          % input -> hidden weights
        wp_kj;          % hidden weights -> output
        w_init_range;   % weight initialization

        % Reinforcement Learning parameters
        beta_p;         % step-size schedule of weight update "learning rate"
        variance;       % variance of perturbation distribution
        varianceRange;
        varDec;
        % covmat;         % action perturbation matrix

        % model state tracking of previous time step
        z_i_prev;       % input layer activation
        z_j_prev;       % hidden layer activation
        z_k_prev;       % output layer activation
        command_prev;   % resulted action

        % model history tracking
        param_num;
        params;
    end

    methods
        function obj = CACLAActor(PARAM)
            obj.input_dim = PARAM{1};
            obj.hidden_dim = PARAM{2};
            obj.output_dim = PARAM{3};
            obj.w_init_range = PARAM{4};
            % obj.wp_ji = rand(obj.output_dim, obj.input_dim) * obj.w_init_range(1); % [0, 1] * w_init_range
            obj.wp_ji = (2 * rand(obj.hidden_dim, obj.input_dim) - 1) * obj.w_init_range(1); % [-1, 1] * w_init_range
            obj.wp_kj = (2 * rand(obj.output_dim, obj.hidden_dim) - 1) * obj.w_init_range(2); % [-1, 1] * w_init_range

            obj.beta_p = PARAM{5};
            obj.varianceRange = PARAM{6};
            obj.variance = obj.varianceRange(1);
            obj.varDec = PARAM{7};
            % obj.covmat = eye(obj.output_dim) * obj.variance;

            obj.param_num = 4;
            obj.params = zeros(1, obj.param_num);

            obj.z_i_prev = zeros(obj.input_dim, 1);
            obj.z_j_prev = zeros(obj.hidden_dim, 1);
            obj.z_k_prev = zeros(obj.output_dim, 1);
            obj.command_prev = 0;
        end

        function update(this)
            % delta_weights(hidden -> output)
            dwp_kj = (this.command_prev - this.z_k_prev) * this.z_j_prev';

            % delta_weights(input -> hidden) [standard backprop]
            dwp_ji = ((1 - this.z_j_prev .^ 2) * this.z_i_prev') * (this.wp_kj * dwp_kj') * this.z_i_prev;

            this.wp_kj = this.wp_kj + this.beta_p * dwp_kj;
            this.wp_ji = this.wp_ji + this.beta_p * dwp_ji * this.z_i_prev';

            % model state tracking
            this.params(1) = sum(sum(abs(this.wp_ji)));
            this.params(2) = sum(sum(this.wp_ji .^ 2));
            this.params(3) = sum(sum(abs(this.wp_kj)));
            this.params(4) = sum(sum(this.wp_kj .^ 2));
        end

        function command = act(this, z_i)
            z_j = tanh(this.wp_ji * z_i);           % activity of hidden layer
            z_k = this.wp_kj * z_j;                 % activity of output layer

            % command = mvnrnd(z_k, this.covmat)';  % perturbation of actor's output multivariate version
            command = mvnrnd(z_k, this.variance)';

            % model state tracking
            this.z_i_prev = z_i;
            this.z_j_prev = z_j;
            this.z_k_prev = z_k;
            this.command_prev = command;
        end

        function command = actHard(this, z_i)
            z_j = tanh(this.wp_ji * z_i);   % activity of hidden layer
            command = this.wp_kj * z_j;     % activity of output layer
        end

        function command = train(this, feature, delta, flag_update)
            if (flag_update && delta > 0)
                this.update();
            end
            command = this.act(feature);
        end
    end
end
