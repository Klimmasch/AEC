classdef CActorG < handle
    properties
        alpha_p; alpha_n;
        input_dim; hidden_dim; output_dim;
        w_init_range;
        wp_ji;      % weights from input to hidden layer
        wp_kj;      % weights from hidden to output layer
        wn_ji;      % weigths for the natural gradient
        cmd_prev;
        z_k_prev;   % previous state of output layer
        z_j_prev;   % hidden ...
        z_i_prev;   % input ...
        covmat;
        type_hidden;
        type_output;
        regulizer;
        variance;
        param_num;
        params;
    end

    methods
        function obj = CActorG(PARAM)
            % PARAM = {Action, alpha_p, alpha_n, lambda, tau, dimension_Feature, initialWeightRange};
            obj.alpha_p = PARAM{1};
            obj.alpha_n = PARAM{2};
            obj.input_dim = PARAM{3};
            obj.hidden_dim = 20;
            obj.output_dim = 1;
            obj.w_init_range = PARAM{4};
            obj.type_hidden = PARAM{5};
            obj.type_output = PARAM{6};
            obj.wn_ji = zeros((obj.input_dim) * obj.hidden_dim + (obj.hidden_dim) * obj.output_dim, 1);
            obj.wp_ji = (2 * rand(obj.hidden_dim, obj.input_dim) - 1) * obj.w_init_range(1);
            obj.wp_kj = (2 * rand(obj.output_dim, obj.hidden_dim) - 1) * obj.w_init_range(2);
            obj.regulizer = 0.005;
            obj.variance = PARAM{7};
            % obj.covmat = diag(ones(obj.output_dim, 1));
            obj.covmat = eye(obj.output_dim) * obj.variance; % indicated no correlation between the output values
            obj.param_num = 7;
            obj.params = zeros(1, obj.param_num);
        end

        function cmd = actDist(this, z_i)
            a_j = this.wp_ji * z_i;         % calculate activity of hidden layer
            % size(a_j) : 20 x 1
            switch this.type_hidden
                case 'tanh'
                    z_j = tanh(a_j);        % normalization from [-Inf,Inf] to [-1,1]
            end
            z_k = this.wp_kj * z_j;         % activity of output layer
            cmd = z_k;
            cmd = mvnrnd(cmd,this.covmat);  % takes z_k as mean and covmat as variance
            cmd = cmd';

            switch this.type_output
                case 'sigmoidal'
                    z_k = 1 ./ (1 + exp(-z_k)); % use sigmoidal activation function %###SWITCH CASE?
                    cmd = 1 ./ (1 + exp(-cmd));
            end

            % save current results
            this.cmd_prev = cmd;
            this.z_i_prev = z_i;
            this.z_j_prev = z_j;
            this.z_k_prev = z_k;
        end

        function z_k = actHard(this, z_i)
            a_j = this.wp_ji * z_i;
            switch this.type_hidden
                case 'tanh'
                    z_j = tanh(a_j);
            end
            z_k = this.wp_kj * z_j; % no addition of noise

            switch this.type_output
                case 'sigmoidal'
                    z_k = 1 ./ (1 + exp(-z_k));
            end
        end

        function update(this, delta)
            delta_k = this.cmd_prev - this.z_k_prev;
            switch this.type_output
                case 'sigmoidal'
                    delta_k = (this.z_k_prev .* (1 - this.z_k_prev)) .* (this.cmd_prev - this.z_k_prev);
            end
            switch this.type_hidden
                case 'tanh'
                    delta_j = (1 - this.z_j_prev .^ 2) .* (this.wp_kj' * delta_k);
            end
            % compatible feature
            % z_j_prev: last activity of hidden layer
            % z_i_prev: last activity of input layer
            dlogp_vp = delta_k * this.z_j_prev';
            dlogp_wp = delta_j * this.z_i_prev';

            % natural gradient
            psi = [dlogp_wp(:); dlogp_vp(:)];
            dwn_ji = delta * psi - (psi' * this.wn_ji) * psi;
            this.wn_ji = this.wn_ji + this.alpha_n * dwn_ji;
            g = this.alpha_p * this.wn_ji;

            % update policy network
            dwp = reshape(g(1 : numel(dlogp_wp)), size(this.wp_ji));            % update for weigths to hidden layer
            dvp = reshape(g(numel(dlogp_wp) + 1 : end), size(this.wp_kj));      % updates for weigths to output layer
            % this.wp_ji = (1 - this.regulizer * this.alpha_p) * this.wp_ji;    % regularization prev: 0.005
            this.wp_ji = this.wp_ji + dwp;
            % this.wp_kj = (1 - this.regulizer * this.alpha_p) * this.wp_kj;    % regularization of output weights (lukas)
            this.wp_kj = this.wp_kj + dvp;

            % save absolute and delta weights
            this.params(1) = sum(sum(abs(this.wp_ji)));
            this.params(2) = sum(sum(this.wp_ji .^ 2));
            this.params(3) = sum(sum(abs(this.wp_kj)));
            this.params(4) = sum(sum(this.wp_kj .^ 2));
            this.params(5) = sum(sum(abs(this.wn_ji)));
            this.params(6) = sum(sum(this.wn_ji .^ 2));
            this.params(7) = norm(g,'fro');

            % this.params(8) = norm(dwp, 'fro');
            % this.params(9) = norm(dvp, 'fro');
            % this.params(10) = psi' * this.wn_ji;
        end

        function cmd = train(this, feature, delta, flag_update)
            if (flag_update)
                this.update(delta);
            end
            cmd = this.actDist(feature);
        end
    end
end
