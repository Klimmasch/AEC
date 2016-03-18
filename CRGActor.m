%%% Continuous Regular Gradient Actor
classdef CRGActor < handle
    properties
        % network parameters
        input_dim;
        output_dim;
        wp_ki;          % input -> hidden weights
        w_init_range;   % weight initialization

        % Reinforcement Learning parameters
        beta_p;         % step-size schedule of weight update "learning rate"
        variance;       % variance of perturbation distribution
        covmat;         % action perturbation matrix

        % model state tracking of previous time step
        z_i_prev;       % input layer activation
        z_k_prev;       % output layer activation
        command_prev;   % resulted action

        % model history tracking
        param_num;
        params;
    end

    methods
        function obj = CRGActor(PARAM)
            obj.input_dim = PARAM{1};
            obj.output_dim = PARAM{2};
            obj.w_init_range = PARAM{3};
            % obj.wp_ki = (2 * rand(obj.output_dim, obj.input_dim) - 1) * obj.w_init_range(1); % [-1, 1] * w_init_range
            obj.wp_ki = rand(obj.output_dim, obj.input_dim) * obj.w_init_range(1); % [0, 1] * w_init_range

            obj.beta_p = PARAM{4};
            obj.variance = PARAM{5};
            obj.covmat = eye(obj.output_dim) * obj.variance;

            obj.param_num = 2;
            obj.params = zeros(1, obj.param_num);

            obj.z_i_prev = zeros(obj.input_dim, 1);
            obj.z_k_prev = zeros(obj.output_dim, 1);
            obj.command_prev = 0;
        end

        function update(this, delta)
            % psi = ((this.command_prev - this.z_k_prev) * this.z_i_prev) / this.variance;
            psi = ((this.z_k_prev - this.command_prev) * this.z_i_prev) / this.variance; % inverse psi
            % histogram(psi,100)
            this.wp_ki = this.wp_ki + this.beta_p * delta * psi';

            % model state tracking
            this.params(1) = sum(sum(abs(this.wp_ki)));
            this.params(2) = sum(sum(this.wp_ki .^ 2));
        end

        function command = act(this, z_i)
            z_k = this.wp_ki * z_i;                 % activity of output layer
            % command = mvnrnd(z_k, this.covmat)';    % perturbation of actor's output multivariate version
            command = mvnrnd(z_k, this.variance)';

            % model state tracking
            this.z_i_prev = z_i;
            this.z_k_prev = z_k;
            this.command_prev = command;
        end

        function command = actHard(this, z_i)
            command = this.wp_ki * z_i;                 % activity of output layer
        end

        function command = train(this, feature, delta, flag_update)
            if (flag_update)
                this.update(delta);
            end
            command = this.act(feature);
        end
    end
end
