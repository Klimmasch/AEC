classdef CActorG < handle
    properties
        alpha_p; alpha_n;
        input_dim; hidden_dim; output_dim;
        w_init_range;
        wp_ji; % weights from input to hidden layer
        wp_kj; % weights from hidden to output layer
        wn_ji; % weigths for the natural gradient
        cmd_prev;
        z_k_prev; % previous state of output layer
        z_j_prev; % hidden ...
        z_i_prev; % input ...
        covmat;
        type_hidden;
        param_num;
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = CActorG(PARAM)
            % PARAM = {Action, alpha_p, alpha_n, lambda, tau, dimension_Feature, initialWeightRange};
            obj.alpha_p = PARAM{1};
            obj.alpha_n = PARAM{2};
            obj.input_dim = PARAM{3};
            obj.hidden_dim = 20;
            obj.output_dim = 2;
            obj.w_init_range = PARAM{4};
            obj.type_hidden = PARAM{5};
            obj.covmat = diag(ones(obj.output_dim, 1)); % indicated no correlation between the output values
            obj.covmat = [.01 0; 0 .01]; % reduce the variance to enable easier learning
            %
            obj.wn_ji = zeros((obj.input_dim)*obj.hidden_dim + (obj.hidden_dim)*obj.output_dim, 1);
            obj.wp_ji = (2*rand(obj.hidden_dim, obj.input_dim)-1) * obj.w_init_range; %*0.15;
            obj.wp_kj = (2*rand(obj.output_dim, obj.hidden_dim)-1) * obj.w_init_range; %*0.22;
            obj.param_num = 6;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function cmd = actDist(this, z_i)
            a_j = this.wp_ji*z_i; % calculate activity of hidden layer
            % size(a_j) : 20 x 1
            switch this.type_hidden
                case 'tanh'
                    z_j = tanh(a_j); % normalization from [-Inf,Inf] to [-1,1]
            end
            z_k = this.wp_kj * z_j; % activity of 1-dimensional output layer
            % now put it into root space (lukas)
%             if (z_k >= 0)
%                 z_k = z_k.^2;
%             else
%                 z_k = -z_k.^2;
%             end

            cmd = mvnrnd(z_k', this.covmat); % takes z_k as mean and covmat as variance
            % cmd = round(cmd, 1); % round to 1 digit (lukas)
            cmd = floor(cmd .* 10) / 10; % round to 1 digit (lukas)
            cmd = cmd';
            %cmd = tanh(cmd)*8 % constrains action space to [-8,8] (lukas)
            % save current results
            this.cmd_prev = cmd;
            this.z_i_prev = z_i;
            this.z_j_prev = z_j;
            this.z_k_prev = z_k;

            cmd = sum(cmd);%make it one-dimensional again for testing purposes (lukas)
        end

        function z_k = actHard(this, z_i)
            a_j = this.wp_ji*z_i;
            switch this.type_hidden
                case 'tanh'
                    z_j = tanh(a_j);
            end
            z_k = this.wp_kj * z_j; %makes it deterministic by excluding covariences/stds/rounding
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function params = update(this, delta)
            params = zeros(1, this.param_num);
            delta_k = this.cmd_prev - this.z_k_prev;
            switch this.type_hidden
                case 'tanh'
                    delta_j = (1-this.z_j_prev.^2).*(this.wp_kj' * delta_k);
            end
            % compatible feature
            % z_j_prev: last activity of hidden layer
            % z_i_prev: last activity of input layer
            dlogp_vp = delta_k * this.z_j_prev';
            dlogp_wp = delta_j * this.z_i_prev';

            % natural gradient
            psi = [dlogp_wp(:);dlogp_vp(:)];
            dwn_ji = delta*psi - (psi'*this.wn_ji)*psi;
            this.wn_ji = this.wn_ji + this.alpha_n * dwn_ji;
            g = this.alpha_p * this.wn_ji;

            % update policy network
            dwp = reshape(g(1:numel(dlogp_wp)), size(this.wp_ji)); % update for weigths to hidden layer
            dvp = reshape(g(numel(dlogp_wp)+1:end), size(this.wp_kj)); % updates for weigths to output layer
            this.wp_ji = (1-0.005*this.alpha_p)*this.wp_ji; % regularization prev: 0.005
            this.wp_ji = this.wp_ji + dwp;
            this.wp_kj = this.wp_kj + dvp;

            % save results
            params(1) = norm(this.wp_ji, 'fro'); params(2) = norm(dwp, 'fro');
            params(3) = norm(this.wp_kj, 'fro'); params(4) = norm(dvp, 'fro');
%             params(3) = norm(this.wp_kj, 'fro');
%             params(4) = norm(this.wp_kj, 'fro');
%             params(5) = psi' * this.wn_ji;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cmd, params] = train(this, feature, delta, flag_update)
            params = zeros(1, this.param_num);
            if (flag_update)
                params = this.update(delta);
            end
            cmd = this.actDist(feature);
        end
    end
end
