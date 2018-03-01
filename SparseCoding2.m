classdef SparseCoding2 < handle
    properties
        nBasis;         % total number of basis
        nBasisUsed;     % number of basis used to encode images in sparse mode
        basisSize;      % size of each (binocular) base vector: patchSize * patchSize * 2 (left + right eye)
        eta;            % learning rate
        temperature;    % temperature in softmax
        basis;          % all basis functions
        basisHist;      % basis functions history
        currentCoef;    % current coeffient matrix
        currentError;   % current reconstruction error
        sizeBatch;      % image batch size's 2nd dimension
    end

    methods
        % Constructor
        % PARAM = [nBasis, nBasisUsed, basisSize, eta, temperature, sizeBatch]
        function obj = SparseCoding2(PARAM)
            obj.nBasis = PARAM(1);
            obj.nBasisUsed = PARAM(2);
            obj.basisSize = PARAM(3);
            obj.eta = PARAM(4);
            obj.temperature = PARAM(5);
            obj.sizeBatch = PARAM(6);

            % obj.basis = rand(obj.basisSize, obj.nBasis) - 0.5;
            % obj.basis = obj.basis * diag(1 ./ sqrt(sum(obj.basis .* obj.basis)));
            % tmpNorm = ones(obj.basisSize, 1) * sqrt(sum(obj.basis .* obj.basis, 1));
            % obj.basis = obj.basis ./ tmpNorm;

            obj.basis = BaseGenerator(0, 0, sqrt(PARAM(3)/2), PARAM(1));
            obj.basisHist = [];

            obj.currentCoef = zeros(obj.nBasis, obj.sizeBatch);     %288x81
            obj.currentError = zeros(obj.basisSize, obj.sizeBatch); %128x81

            %TODO maybe reimplement reloading basis functions
        end

        %%% Encode the input images accoring to softmax distribution
        %   @param imageBatch:  input image patches batch
        function softmaxEncode(this, imageBatch)
            this.currentCoef = this.currentCoef * 0; % needs to be tested if resetting is necessary
            tmp = imageBatch;
            for count = 1 : this.nBasisUsed
                corrl = abs(this.basis' * tmp) / this.temperature;
                corrl = corrl - kron(ones(this.nBasis, 1), max(corrl));
                softmaxcorr = softmax(corrl);

                softmaxcorr = tril(ones(this.nBasis)) * softmaxcorr - kron(ones(this.nBasis, 1), rand(1, this.sizeBatch));
                softmaxcorr(softmaxcorr < 0) = 2;
                [~, index] = min(softmaxcorr);
                corrl = this.basis' * tmp;
                linearIndex = sub2ind(size(corrl), index, 1 : this.sizeBatch);
                this.currentCoef(linearIndex) = this.currentCoef(linearIndex) + corrl(linearIndex);
                tmp = imageBatch - this.basis * this.currentCoef;
            end
            this.currentError = tmp;
        end

        %%% Encode the input images with the best matched basis
        %   @param imageBatch:  input image patches batch
        function sparseEncode(this, imageBatch)
            % batch_size = size(Images, 2);
            % Coef = zeros(this.nBasis, batch_size);
            % I = Images;
            % for (count = 1:this.nBasisUsed)
            %     corrl = this.Basis'*I;
            %     [~, index] = max(abs(corrl));
            %     alpha = diag(this.Basis(:, index)'*I);
            %     linearIndex = sub2ind(size(corrl), index, 1:batch_size);
            %     Coef(linearIndex) = Coef(linearIndex) + alpha';
            %     I = Images - this.Basis*Coef;
            % end
            % Error = I;

            this.currentCoef = this.currentCoef * 0;                                    % needs to be tested if resetting is necessary
            corrl = this.basis' * imageBatch;                                           % correlation of each basis with each patch
            corrBB = this.basis' * this.basis;                                          % correlation between basis
            for count = 1 : this.nBasisUsed
                [~, index] = max(abs(corrl));                                           % indices of bases with max correlation per patch
                linearIndex = sub2ind(size(corrl), index, 1 : this.sizeBatch);          % corresponding linear indicies in correlation matrix
                pCorr = corrl(linearIndex);                                             % vector of correlations per patch (coefs per patch)
                this.currentCoef(linearIndex) = this.currentCoef(linearIndex) + pCorr;  % calculate new correlation coefficients
                corrl = corrl - bsxfun(@times, corrBB(:, index), pCorr);                % (see Yu's doc)
            end
            this.currentError = imageBatch - this.basis * this.currentCoef;             % 128x81 = 128x81 - 128x288 * 288x81
        end

        %%% Calculate the correlation between input image and the basis
        %   @param imageBatch:  input image patches batch
        function fullEncode(this, imageBatch)
            this.currentCoef = this.basis' * imageBatch;
            this.currentError = imageBatch - this.basis * this.currentCoef;
        end

        %%% Update basis functions
        function stepTrain(this)
            deltaBases = this.currentError * this.currentCoef' / size(this.currentError, 2);
            this.basis = this.basis + this.eta * deltaBases;
            this.basis = bsxfun(@rdivide, this.basis, sqrt(sum(this.basis .^ 2)));
        end

        %%% Track the evolution of all basis functions over time
        function saveBasis(this)
            this.basisHist = cat(3, this.basisHist, this.basis);
        end

        %%% Display the Basis functions (Zhao Yu code) at iteration t
        function displayBasis(this, t)
            %how to arrange the basis (rows, col)
            R = 16;
            C = 18;
            len = 1;
            % basisTrack = this.drecord.basisTrack(1:len);
            basisTrack{1} = this.basis;
            %checkPoint = 1;

            endBasis = basisTrack{end}(1 : end / 2, :);
            leftEnergy = abs(sum(endBasis .^ 2) - 0.5);
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
            for j = 1 : len
                A = basisTrack{j}(:, I);
                % B = reshape(A, di*sqrt(num/2), sqrt(num/2)*2);
                B = reshape(A, di * R, C);
                B = B / max(max(abs(B))) + 0.5;
                C = padarray(padarray(blockproc(B, [di, 1], fun1) - 1, [1 1],'post') + 1,[2, 2]);
                imshow(C);
                % title(num2str(checkPoint(j)));
                title(num2str(t));
                drawnow;
            end
        end
    end
end
