classdef SparseCoding < handle
    properties
        Basis_num_used; %number of basis used to encode in sparse mode
        Basis_size;     %size of each base vector
        Basis_num;      %total basis number
        Basis;          %all the basis
        Basis_hist;     %save basis at regular intervals
        eta;            %learning rate
        Temperature;    %temperature in softmax
        Dsratio;        %Downsampling ratio (to produce 8x8)
    end

    methods
        %PARAM = {Basis_num_used, Basis_size, Basis_num, eta, Temperature, Dsratio, Basis, loadBasis};
        function obj = SparseCoding(PARAM)
            obj.Basis_num_used = PARAM{1};
            obj.Basis_size = PARAM{2};
            obj.Basis_num = PARAM{3};
            obj.eta = PARAM{4};
            obj.Temperature = PARAM{5};
            obj.Dsratio = PARAM{6};
            if (PARAM{8})
                obj.Basis = PARAM{7};
                obj.Basis_hist = PARAM{7};
            else
                a = rand(obj.Basis_size, obj.Basis_num)-0.5; % basis function set
                a = a*diag(1./sqrt(sum(a.*a)));
                thenorm = ones(obj.Basis_size, 1)*sqrt(sum(a.*a, 1));
                a = a./thenorm;
                obj.Basis = a;
                obj.Basis_hist = a;
            end
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
        function [Coef, Error] = softmaxEncode(this, Images)
            batch_size = size(Images, 2);
            Coef = zeros(this.Basis_num, batch_size);
            I = Images;
            for count = 1:this.Basis_num_used
                corr = abs(this.Basis'*I)/this.Temperature;
                corr = corr - kron(ones(this.Basis_num, 1), max(corr));
                softmaxcorr = softmax(corr);

                softmaxcorr = tril(ones(this.Basis_num)) * softmaxcorr - kron(ones(this.Basis_num, 1), rand(1, batch_size));
                softmaxcorr(softmaxcorr<0) = 2;
                [~, index] = min(softmaxcorr);
                corr = this.Basis'*I;
                linearindex = sub2ind(size(corr), index, 1:batch_size);
                Coef(linearindex) = Coef(linearindex) + corr(linearindex);
                I = Images - this.Basis*Coef;
            end
            Error = I;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Encode the input images with the best matched basis
        %%%
        %%% Images are the input images batch
        %%%
        %%% Coef is the output Coefficients
        %%% Error is the reconstructin error using current coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [coef, error] = sparseEncode(this, imageBatch)
            % batch_size = size(Images, 2);
            % Coef = zeros(this.Basis_num, batch_size);
            % I = Images;
            % for (count = 1:this.Basis_num_used)
            %     corr = this.Basis'*I;
            %     [~, index] = max(abs(corr));
            %     alpha = diag(this.Basis(:, index)'*I);
            %     linearindex = sub2ind(size(corr), index, 1:batch_size);
            %     Coef(linearindex) = Coef(linearindex) + alpha';
            %     I = Images - this.Basis*Coef;
            % end
            % Error = I;
            size_Batch = size(imageBatch, 2);
            coef = zeros(this.Basis_num, size_Batch);
            corr = this.Basis'*imageBatch;      %correlation of each basis with each patch
            corrBB = this.Basis'*this.Basis;    %correlation between basis
            for count = 1:this.Basis_num_used
                [~, index] = max(abs(corr));                            %indices of bases with max corr per patch
                linearindex = sub2ind(size(corr), index, 1:size_Batch); %corresponding linear indices in corr matrix
                pCorr = corr(linearindex);                              %vector of correlations per patch (coefs per patch)
                coef(linearindex) = coef(linearindex) + pCorr;          %stores corr coefs into coef matrix
                corr = corr - bsxfun(@times, corrBB(:, index), pCorr);  %(see Yu's doc)
            end
            error = imageBatch - this.Basis*coef;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate the correlation between input image and the basis
        %%%
        %%% Images are the input image batch
        %%%
        %%% Coef is the output correlation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Coef, Error] = fullEncode(this, Images)
            Coef = this.Basis'*Images;
            Error = Images - this.Basis*Coef;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Update the basis
        %%%
        %%% Coef is the input coefficient
        %%% Error is the input error
        %%% debugmode indicates whether some intermedia should be recorded;
        %%%
        %%% Basis_Change is the changing amount of the basis in current
        %%% update
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function updateBasis(this, coef, error)
            % da = Error * Coef';
            % da = da/size(Error, 2);
            % basis = this.Basis + this.eta*da;
            % basis = basis./kron(ones(this.Basis_size, 1), sqrt(sum(basis.^2)));
            % this.Basis = basis;

            deltaBases = error * coef' / size(error, 2);
            this.Basis = this.Basis + this.eta * deltaBases;
            this.Basis = bsxfun(@rdivide, this.Basis, sqrt(sum(this.Basis.^2)));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% train sparse coding for one step
        %%%
        %%% Images is the input image batch
        %%% debugmode indicates whether some intermedia should be recorded;
        %%%
        %%% Error is the reconstruction error using the best matched coefficients
        %%% Basis_picked indicates which basis are picked to encode
        %%% Basis_Entropy is the entropy of each base
        %%% Basis_Change is the changing amount of the basis in current
        %%% update
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [error, coef] = stepTrain(this, Images)
            % [Coef, Error] = this.softmaxEncode(Images);
            [coef, error] = this.sparseEncode(Images); %Matching Pursuit
            updateBasis(this, coef, error); %Gradient descent (Basis change)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the parameters in a file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveClass(this, configfile)
            Basis = this.Basis;
            save(configfile,'Basis','-append');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% save the Basis during training
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function saveBasis(this)
            this.Basis_hist = cat(3, this.Basis_hist, this.Basis);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Display the Basis functions (Zhao Yu code) at iteration t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function displayBasis(this, t)
            %how to arrange the basis (rows, col)
            R = 16;
            C = 18;
            len = 1;
            % basisTrack = this.drecord.basisTrack(1:len);
            basisTrack{1} = this.Basis;
            %checkPoint = 1;

            endBasis = basisTrack{end}(1:end/2,:);
            leftEnergy = abs(sum(endBasis.^2)-0.5);
            [~, I] = sort(leftEnergy);

            % h = gcf;
            % set(h,'Position',[1 1 800 600]);
            % scrsz = get(0,'ScreenSize');
            % set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
            subplot(1, 1, 1);
            [di, num] = size(basisTrack{1});
            fun1 = @(blc_struct) padarray(padarray(reshape(permute(padarray(reshape(blc_struct.data, sqrt(di / 2), ...
                     sqrt(di / 2), 2), [1, 1], 'pre'), [1, 3, 2]), (sqrt(di / 2) + 1) * 2, sqrt(di / 2) + 1), ...
                     [1, 1], 'post') - 1, [1 1], 'pre') + 1;
            for j = 1:len
                A = basisTrack{j}(:, I);
                % B = reshape(A, di*sqrt(num/2), sqrt(num/2)*2);
                B = reshape(A, di*R, C);
                B = B/max(max(abs(B))) + 0.5;
                C = padarray(padarray(blockproc(B,[di, 1], fun1)-1,[1 1],'post')+1,[2, 2]);
                imshow(C);
                % title(num2str(checkPoint(j)));
                title(num2str(t));
                drawnow;
            end
        end
    end
end
