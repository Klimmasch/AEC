classdef CCriticG < handle
    properties
        alpha_v;
        eta;
        gamma;
        input_dim;
        v_ji; %weights
        v_init_range;
        z_i_prev;
        z_j_prev;
        J;
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% initialization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = CCriticG(PARAM)
            % PARAM = {alpha_v, eta, gamma, input_dim, v_init_range};
            obj.alpha_v = PARAM{1};
            obj.eta = PARAM{2}; % compute J
            obj.gamma = PARAM{3}; % discount
            obj.input_dim = PARAM{4}+1;%###!
            obj.v_init_range = PARAM{5};
            %
            obj.v_ji = (2*rand(1, obj.input_dim)-1)*obj.v_init_range;%0.12;
            obj.z_j_prev = 0;
            obj.J = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calculates the output of the network from the input feature
        %%% vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z_j = forward(this, z_i)
            z_j = this.v_ji * z_i;%[z_i;1];%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calculates gradient and updates the weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [delta, params] = update(this, z_j, reward)
            params = zeros(1, 2);
            % update network
            this.J = (1-this.eta) * this.J + this.eta*reward;%0;
            delta = reward - this.J + this.gamma * z_j - this.z_j_prev;
            dv_ji = this.alpha_v * delta * this.z_i_prev';
            this.v_ji = this.v_ji +  dv_ji;
            params(1) = norm(this.v_ji, 'fro'); %2 norm of column vector
            params(2) = delta;%z_i;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% inputs: feature vector, reward
        %%% outputs: gradient for the actor (delta), [normalized weights of
        %%% network, delta again]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [delta, params] = train(this, z_i, reward, flag_update)
            delta = 0;
            params = [0, 0];
            z_j = this.forward(z_i); % calculate internal activity (feature times weights)
            if (flag_update)
                [delta, params] = this.update(z_j, reward);
            end
            this.z_i_prev = z_i;%[z_i;1];%
            this.z_j_prev = z_j;

            %display(sprintf('critic: \n\t z_j: %03f,\n\t z_j_prev: %03f,\n\t z_i: %d,\n\t z_i_prev: %d,\n\t reward:%f, \n\t delta:%03f',z_j, this.z_j_prev, mean(z_i), mean(this.z_i_prev),reward, delta))
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
