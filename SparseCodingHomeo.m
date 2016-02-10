classdef SparseCodingHomeo < handle
    properties        
        Basis;%all the basis
        Basis_hist; %save basis at regular intervals
       
        Basis_num_used;%number of basis used to encode in sparse mode
        Basis_size;%size of each base vector
        Basis_num;%total basis number       
        eta;%learning rate 
        Temperature;%temperature in softmax
        Dsratio;    %Downsampling ratio (to produce 8x8)
        
        switch_sym; % we use ON/OFF symmetry
        batch_size; % how many 'imagelets' do we take for each learning step?
        display_every; % delay between snapshots (to the screen or movie)
        
        % learning parameters
        frac; % we take *at most* frac active filters in the Matching Pursuit during the lerning phase
        noise_var_ssc;  % relative threshold for SSC corresponding to an
                           % estimate of the ratio of background
                           % noise over total signal energy
        var_eta_ssc;    % used to ensure that all filters
        switch_choice;
        n_quant;
        switch_Mod;

        where;
        
        
        %Variables
        gain_rand;
        Pz_j;
        Mod;
        Pz_j_;
        
        dA;
    end
    methods
        %PARAM = {Basis_num_used,Basis_size,Basis_num,eta,Temperature,Dsratio,Basis_S,loadBasis};
        function obj = SparseCodingHomeo(PARAM)
            obj.Basis_num_used = PARAM{1};
            obj.Basis_size = PARAM{2};
            obj.Basis_num = PARAM{3};
            obj.eta = PARAM{4};
            obj.Temperature = PARAM{5};
            obj.Dsratio = PARAM{6};
            if(PARAM{8})
                obj.Basis = PARAM{7};
                obj.Basis_hist = PARAM{7};
                
            else
                a=rand(obj.Basis_size,obj.Basis_num)-0.5; % basis function set
                a=a*diag(1./sqrt(sum(a.*a)));
                obj.gain_rand = sqrt(sum(a.*a))';
                thenorm = ones(obj.Basis_size,1)*sqrt(sum(a.*a,1));
                a=a./thenorm;
                obj.Basis = a;
                obj.Basis_hist = a;
            end
            obj.switch_sym = 1;
            obj.batch_size = 100;
            obj.display_every = 50;
            obj.frac = 0.25;
            obj.noise_var_ssc = 0.002;
            obj.var_eta_ssc = 1/20;
            obj.switch_choice = obj.var_eta_ssc > 0;
            obj.n_quant = 512;%???
            obj.switch_Mod = 1;
            
            obj.Pz_j=1/obj.n_quant*ones(obj.n_quant,obj.Basis_num); 
            obj.Mod=cumsum(obj.Pz_j);
            
            obj.where = ['../results/' datestr(now, 30)];            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% encode the image accoring to softmax distribution
        %%%
        %%% Images is the batch input
        %%% debugmode indicates whether some intermedia should be recorded;
        %%%
        %%% Coef is the output Coefficients for each basis and images
        %%% Error is the reconstruction error using current coefficients
        %%% Basis_picked indicates which basis are picked to encode
        %%% Basis_Entropy is the entropy of each base
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Coef, Error] = softmaxEncode(this,Images)
            batch_size = size(Images,2);
            Coef = zeros(this.Basis_num,batch_size);
            I = Images;
            for count = 1:this.Basis_num_used
                corr = abs(this.Basis'*I)/this.Temperature;
                corr = corr - kron(ones(this.Basis_num,1),max(corr));
                softmaxcorr = softmax(corr);
                
                softmaxcorr = tril(ones(this.Basis_num))*softmaxcorr - kron(ones(this.Basis_num,1),rand(1,batch_size));
                softmaxcorr(softmaxcorr<0) = 2;
                [~,index] = min(softmaxcorr);
                corr = this.Basis'*I;
                linearindex = sub2ind(size(corr),index,1:batch_size);
                Coef(linearindex) = Coef(linearindex) + corr(linearindex);
                I = Images - this.Basis*Coef;
            end
            Error = I;
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Encode the input images with the best matched basis
        %%% 
        %%% Images are the input images batch
        %%% A_rand = rand(e.L,e.M)-0.5; A_rand = A_rand*diag(1./sqrt(sum(A_rand.*A_rand)));

        %%% Coef is the output Coefficients
        %%% Error is the reconstructin error using current coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [coef,error] = sparseEncode(this,imageBatch)
            %             batch_size = size(Images,2);
            %             Coef = zeros(this.Basis_num,batch_size);
            %             I = Images;
            %             for count = 1:this.Basis_num_used
            %                 corr = this.Basis'*I;
            %                 [~,index] = max(abs(corr));
            %                 alpha = diag(this.Basis(:,index)'*I);
            %                 linearindex = sub2ind(size(corr),index,1:batch_size);
            %                 Coef(linearindex) = Coef(linearindex) + alpha';
            %                 I = Images - this.Basis*Coef;
            %             end
            %             Error = I;
            size_Batch = size(imageBatch,2); %X = imageBatch
            this.batch_size = size_Batch;
            coef = zeros(this.Basis_num,size_Batch); % initialize coeffs for LGM
            dA = zeros(size(this.Basis)); % initalize weight's gradient
            
            % --------------------------------------------
            % SPARSIFICATION
            [coef, dA, Pz_j_]=mp_fitS(this.Basis,imageBatch,this.noise_var_ssc,this.frac,this.switch_choice,this.Mod,0,0,this.switch_sym);%
            this.Pz_j_ = Pz_j_;
            this.dA = dA;
            % residual
            error = imageBatch - this.Basis*coef;
            
        end
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate the corre lation between input image and the basis
        %%% 
        %%% Images are the input image batch
        %%% 
        %%% Coef is the output correlation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Coef,Error] = fullEncode(this,Images)
            Coef = this.Basis'*Images;
            Error = Images - this.Basis*Coef;
        end              
        
        
        function updateBasis(this,coef,error)
            %% --------------------------------------------
            %% apply the learning gradient (dA) on the network (n) and modify
            %% homeostatic variables

                eta = this.eta;% sparsify by Matching Pursuit
            
            %% LEARNING : it is the same for both methods
            if (eta>0), %        
                % 1) updates basis functions

                %% if you increase the batch size, the gradient increases proportionnally
                %% (with ergodicity...), so by Knuth programming law ("Thou shall
                %% make your program scale invariant")
                dA = this.dA;
                dA = dA/this.batch_size;

                %% applies the gradient descent
                this.Basis = this.Basis + eta * dA;%

            end% end learning loop


            %% HOMEOSTASIS : it's different for both methods
            % 2) update the norm and average use of every neuron (homeostatic rules)
            normA=sqrt(sum(this.Basis.*this.Basis));
            % normalization
            for i_M=1:this.Basis_num %over basis functions
                this.Basis(:,i_M)=this.Basis(:,i_M)/normA(i_M);
            end

            if this.var_eta_ssc >0,
                %  adaptive rule for homeostasis
                t_homeo = 1/this.var_eta_ssc; % TODO: remove? min(t,1/e.var_eta_ssc); %
                this.Pz_j=(1-1/t_homeo)*this.Pz_j+1/t_homeo*this.Pz_j_;% update statistics
                this.Mod=cumsum(this.Pz_j); %
            end        
            
        
        end
        
        function [error,coef] = stepTrain(this,Images)
            %             [Coef,Error] = this.softmaxEncode(Images);
            [coef,error] = this.sparseEncode(Images);    %Matching Pursuit
            updateBasis(this,coef,error);           %Gradient descent (Basis change)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the parameters in a file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveClass(this,configfile)
            Basis = this.Basis;
            save(configfile,'Basis','-append');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the Basis during training
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveBasis(this)
            this.Basis_hist = cat(3,this.Basis_hist,this.Basis);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Display the Basis functions (Zhao Yu code) at iteration t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function displayBasis(this,t)
            R = 16; C = 18;      %how to arrange the basis (rows, col)
            len = 1;
            % basisTrack = this.drecord.basisTrack(1:len);
            basisTrack{1} = this.Basis;
            %checkPoint = 1;
            
            endBasis = basisTrack{end}(1:end/2,:);
            leftEnergy = abs(sum(endBasis.^2)-0.5);
            [~,I] = sort(leftEnergy);
            
%             h = gcf;
%             set(h,'Position',[1 1 800 600]);
%             scrsz = get(0,'ScreenSize');
%             set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
            subplot(1,1,1);
            [di,num] = size(basisTrack{1});
            fun1 = @(blc_struct) padarray(padarray(reshape(permute(padarray(reshape(blc_struct.data,sqrt(di/2),sqrt(di/2),2),[1,1],'pre'),[1,3,2]),(sqrt(di/2)+1)*2,sqrt(di/2)+1),[1,1],'post')-1,[1 1],'pre')+1;
            for j = 1:len
                A = basisTrack{j}(:,I);
                %                 B = reshape(A,di*sqrt(num/2),sqrt(num/2)*2);
                B = reshape(A,di*R,C);
                B = B/max(max(abs(B))) + 0.5;
                C = padarray(padarray(blockproc(B,[di,1],fun1)-1,[1 1],'post')+1,[2,2]);
                imshow(C);
%                 title(num2str(checkPoint(j)));
                title(num2str(t));
                
                drawnow;
            end
        end
    end 
end
