%%%
% Continuous Actor Critic Learning Automaton with variance Actor
% Lukas's interpretation of CACLA appoach
%%%
classdef CACLAVarActorLu < handle
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
        covmat;         % action perturbation matrix

        % model state tracking of previous time step
        z_i_prev;       % input layer activation
        z_j_prev;       % hidden layer activation
        z_k_prev;       % output layer activation
        command_prev;   % resulted action
        deltaVar;       % variance of delta signal
        eta;            % scaling factor of delta variance
        updateCount;

        % model history tracking
        param_num;
        params;
        
        regularizer;
    end

    methods
        function obj = CACLAVarActorLu(PARAM)
            obj.input_dim = PARAM{1}(1);
            obj.hidden_dim = PARAM{1}(2);
            obj.output_dim = PARAM{1}(3);
            obj.w_init_range = PARAM{2};
            % obj.wp_ji = rand(obj.output_dim, obj.input_dim) * obj.w_init_range(1); % [0, 1] * w_init_range
            obj.wp_ji = (2 * rand(obj.hidden_dim, obj.input_dim) - 1) * obj.w_init_range(1); % [-1, 1] * w_init_range
            obj.wp_kj = (2 * rand(obj.output_dim, obj.hidden_dim) - 1) * obj.w_init_range(2); % [-1, 1] * w_init_range

            obj.beta_p = PARAM{3};
            obj.varianceRange = PARAM{4};
            obj.variance = obj.varianceRange(1);
            obj.deltaVar = PARAM{5};
            obj.eta = PARAM{6};
            obj.varDec = PARAM{7};
            obj.covmat = eye(obj.output_dim) * obj.variance;

            obj.param_num = 5;
            obj.params = zeros(1, obj.param_num);

            obj.z_i_prev = zeros(obj.input_dim, 1);
            obj.z_j_prev = zeros(obj.hidden_dim, 1);
            obj.z_k_prev = zeros(obj.output_dim, 1);
            obj.command_prev = zeros(obj.output_dim, 1);
            obj.updateCount = 0;
            obj.regularizer = 0.999; % percentage of the weights the vector is scaled to
        end

        function update(this, delta)
            this.updateCount = ceil(delta / sqrt(this.deltaVar));
            
            % delta_weights(hidden -> output)
            dwp_kj = (this.command_prev - this.z_k_prev) * this.z_j_prev';

            % delta_weights(input -> hidden)
            tmpVector = ((this.command_prev - this.z_k_prev)' * this.wp_kj)';
            dwp_ji = ((1 - this.z_j_prev .^ 2) * this.z_i_prev') .* repmat(tmpVector, 1, this.input_dim);
            
%             this.prams(5) = (this.beta_p * dwp_kj) * this.updateCount;
%             this.wp_kj = this.wp_kj + this.params(5);

            this.wp_kj = this.wp_kj + (this.beta_p * dwp_kj) * this.updateCount;
            
            this.params(4) = mean(mean(abs((this.beta_p * dwp_ji) * this.updateCount)));  % this tracks the change in weight because of updates
            this.params(5) = mean(mean(abs((this.regularizer * this.wp_ji) - this.wp_ji))); % this tracks changes due to regularization
            
            this.wp_ji = (this.regularizer * this.wp_ji) + (this.beta_p * dwp_ji) * this.updateCount;
        end

        function command = act(this, z_i)
            z_j = tanh(this.wp_ji * z_i);           % activity of hidden layer
            z_k = this.wp_kj * z_j;                 % activity of output layer

            command = mvnrnd(z_k, this.covmat)';    % perturbation of actor's output by multivariate Gaussian

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
            % approximate variance in TD signal
            this.deltaVar = (1 - this.eta) * this.deltaVar + this.eta * delta ^ 2;
            if (flag_update && delta > 0)
                this.update(delta);
            else
                this.updateCount = 0;
            end
            command = this.act(feature);

            % model state change tracking
            this.params(1) = sum(sum(abs(this.wp_ji)));
            this.params(2) = sum(sum(abs(this.wp_kj)));
            this.params(3) = this.updateCount;
        end
    end
end
