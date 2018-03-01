% clear;


%Generates set of bases with specific frequency and orientation distribution
%% ####################################################################################################################   

function basis = BaseGenerator(saveBasis, displayBasis, baseSize, number)
    % s=8; %base size
    % n=300; %number of bases to generate

    s = baseSize;
    n = number;
    
    sigma=2.5; %sigma of gaussian envelop of gabor patch
    freq_sig=1; %sigma of exponential for wavelength distribution of bases

    % a=rand(s*s*2,n)-0.5; %initialize bases with random weights and zero mean
    % a=a*diag(1./sqrt(sum(a.*a))); %normalize to unit variance
    % thenorm = ones(obj.Basis_size,1)*sqrt(sum(a.*a,1));
    % a=a./thenorm;


    a=[];
    for i=1:n
        orient=rand()*180; %orientation meassured from horizontal axis anticlockwise in degrees
        theta=rand()*360; %phase shift of sin wave in degrees
        % lambda=rand()*s*3/3+s*3/3; % wavelength of sin wave in px
        % lambda = 1/2*s;
        lambda=rand()*s*5/3+s*1/3;

        shift_x=rand()*s*2/3-s*1/3; %shift wavelet away from patch center
        shift_y=rand()*s*2/3-s*1/3;

    %     data=laprnd(1000000, 1, 0, 0.5);
    %     data=data(abs(data)<1);
    %     hist(data,-1:0.05:1);
        delta_orient=2;
        while(abs(delta_orient)>=1)
            delta_orient=laprnd(1, 1, 0, 1);
        end
        delta_orient=delta_orient*90;

        dispa=rand()*s/4-s/8;
    %     dispa=0;

        base_l = Gabor(sigma,orient,lambda,theta,lambda/(0.8*s),s,[shift_x shift_y],+dispa);
        base_r = Gabor(sigma,orient+delta_orient,lambda,theta,lambda/(0.8*s),s,[shift_x shift_y],-dispa);

        %normalize left and right part of base
        base_l=base_l(:)-mean(base_l(:));
        base_r=base_r(:)-mean(base_r(:));
        base_l=base_l/std(base_l); %normalize to unit variance
        base_r=base_r/std(base_r); %normalize to unit variance

        base=[base_l(:);base_r(:)];
        a=[a, base];
    end

    basis=a;

    basis=bsxfun(@rdivide,basis,sqrt(sum(basis.^2)));

    if saveBasis
        save('basis.mat','basis');
    end


    %% ####################################################################################################################   
    %Display Bases
    
    if displayBasis
        
        %####################################
        orange=[255 153 24]./255;
        blue=[39 48 73]./255;
        green=[146 208 80]./255;
        lightblue=[70 166 199]./255;
        purple=[139 84 166]./255;
        brown=[167 121 29]./255;
        gray=[80 80 80]/255;
        %####################################
    
        % load('model_2.mat');
        % basis=model.scmodel.Basis;

        %select specific subset of basis
        subset = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
        % subset = [1 2 3 4 5 6];
        % bgcolor= orange;
        % bgcolor= [195 117 18]./255;
        bgcolor= [255 255 255]./255;

        basisDisp = basis(:,subset);

        s=sqrt(size(basisDisp,1)/2);
        n=size(basisDisp,2);

        block=mat2gray(col2im(basisDisp,[s 2*s],[n*s 2*s],'distinct'));

        M=3;%5; %number of bases displayed in one row
        N=2;%3; %maximum number of rows to plot

        k=1;
        pic=[];

        for j=1:N
            row=[];
            for i=1:M
                base_l=block(1+(k-1)*s:k*s,1:s);
                base_r=block(1+(k-1)*s:k*s,s+1:2*s);
                brick=[-ones(1,s); base_l; -ones(1,s); base_r; -ones(1,s)]; %left part on top
                brick=[-ones(2*s+3,1) brick -ones(2*s+3,1)];
                row=[row brick];
                k=k+1;
            end

            pic=[pic;row];
        end

        % %add frame
        % pic=[-ones(size(pic,1),1) pic -ones(size(pic,1),1)];
        % pic=[-ones(1,size(pic,2)); pic; -ones(1,size(pic,2))];

        imshow(pic);

        %colorize
        red = pic;
        green = pic;
        blue = pic;

        red(red==-1)= bgcolor(1);
        green(green==-1)= bgcolor(2);
        blue(blue==-1)= bgcolor(3);

        test = zeros(size(pic,1),size(pic,2),3);

        test(:,:,1)=red;
        test(:,:,2)=green;
        test(:,:,3)=blue;

        imshow(test)
    end
end

function y  = laprnd(m, n, mu, sigma)
    %LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
    %   with mean mu and standard deviation sigma. 
    %   mu      : mean
    %   sigma   : standard deviation
    %   [m, n]  : the dimension of y.
    %   Default mu = 0, sigma = 1. 
    %   For more information, refer to
    %   http://en.wikipedia.org./wiki/Laplace_distribution

    %   Author  : Elvis Chen (bee33@sjtu.edu.cn)
    %   Date    : 01/19/07

    %Check inputs
    if nargin < 2
        error('At least two inputs are required');
    end

    if nargin == 2
        mu = 0; sigma = 1;
    end

    if nargin == 3
        sigma = 1;
    end

    % Generate Laplacian noise
    u = rand(m, n)-0.5;
    b = sigma / sqrt(2);
    y = mu - b * sign(u).* log(1- 2* abs(u));
end



function gb=Gabor(sigma,orient,wavel,phase,aspect,pxsize,shift,dispa)
    %function gb=gabor(sigma,orient,wavel,phase,aspect)
    %
    % This function produces a numerical approximation to 2D Gabor function.
    % Parameters:
    % sigma  = standard deviation of Gaussian envelope, this in-turn controls the
    %          size of the result (pixels)
    % orient = orientation of the Gabor clockwise from the vertical (degrees)
    % wavel  = the wavelength of the sin wave (pixels)
    % phase  = the phase of the sin wave (degrees)
    % aspect = aspect ratio of Gaussian envelope (0 = no "width" to envelope, 
    %          1 = circular symmetric envelope)
    % pxsize = the size of the filter (optional). If not specified, size is 5*sigma.

    if nargin<6
      pxsize=fix(5*sigma);
    end

    % if mod(pxsize,2)==0, pxsize=pxsize+1; end %give mask an odd dimension
    if mod(pxsize,2)~=0
        [x, y]=meshgrid(-fix(pxsize/2):fix(pxsize/2),-fix(pxsize/2):fix(pxsize/2));
    else
        [x, y]=meshgrid(-fix(pxsize/2)+0.5:fix(pxsize/2)-0.5,-fix(pxsize/2)+0.5:fix(pxsize/2)-0.5);
    end

    % Add Disparity
    x = x+dispa;

    % Rotation 
    orient=-orient*pi/180;
    x_theta=x*cos(orient)+y*sin(orient);
    y_theta=-x*sin(orient)+y*cos(orient);

    phase=phase*pi/180;
    freq=2*pi./wavel;

    shift_x=shift(1);
    shift_y=shift(2);

    gb=exp(-0.5.*( (((x_theta+shift_x).^2)/(sigma^2)) ...
                 + (((y_theta+shift_y).^2)/((aspect*sigma)^2)) )) ...
       .* (cos(freq*y_theta+phase) - cos(phase).*exp(-0.25.*((sigma*freq).^2)));

    %% Print error if wavelet does not sum up to zero

    % if abs(sum(sum(gb)))./sum(sum(abs(gb)))>0.01
    %   disp('WARNING unbalanced Gabor');
    %   abs(sum(sum(gb)))./sum(sum(abs(gb)))
    % end
end