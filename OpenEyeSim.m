function OpenEyeSim( train_time, randomization_i )
%Main script for launching application
%   learned_file, file with policy and sparse coding to test if not empty
%   texture_file: texture settings files
%   development: 0 - normal, 1 - 8 months, 2 - from 1 to 8 months
%   development_time in iterations (freeze policy)
%   train_time in number of iterations
%   sparse_coding

learned_file = '';texture_file = 'Textures_New.mat';sparse_coding = 'nonhomeo';

rng(randomization_i);

%Init model
model = config(learned_file, texture_file, train_time, sparse_coding); %instantiate model object
model_data = Model_test_data();   %create instance of object

%Predefine testing phases
testing_window = 0;%1000;
testing_counter = 0;

%Testing flag
phase = 'testing' % testing or training, we do this every 10% of learning

savepath = sprintf('model_%s_%i_%i_%i_%s_%i.mat',datestr(now),train_time,sparse_coding,randomization_i);
savepath_tests = sprintf('model_tests_%s_%i_%i_%i_%s_%i.mat',datestr(now),train_time,sparse_coding,randomization_i);

%Control Parameters
vergemax = 16;
anglenew = 0;     %init new vergence command
tmpX = zeros(3);

t = 0;  % # of iterations
status_done = 0;

%image process vars
patch_size = 8;

Dsratio_L = model.scmodel_Large.Dsratio;     %downsampling ratio (Large scale)
Dsratio_S = model.scmodel_Small.Dsratio;     %downsampling ratio (Small scale)

fovea = [128 128];
fovea_L = patch_size + patch_size^2/2^log2(Dsratio_L);    %fovea size (Large scale)
fovea_S = patch_size + patch_size^2/2^log2(Dsratio_S);    %fovea size (Small scale)

s_L = patch_size/Dsratio_L;         %steps of overlap in the ds image
s_S = patch_size/Dsratio_S;         %steps of overlap in the ds image

nc_L = fovea_L-patch_size+1;                %number of patches per column (slide of 1 px)
nc_S = fovea_S-patch_size+1;                %number of patches per column (slide of 1 px)

h1 = 240; w1 = 320; pixSize = 1;                 %input image size

%Camera parameters
offset = 0;    %vertical offset between left and right (0 in the Simulator!!!)
f = 257.34;    %focal length [px]
baseline = 0.056;  %interocular distance (baseline)

%Textures
texture_path = sprintf('config/%s',texture_file);
texture = load(texture_path);texture = texture.texture;
n_textures = length(texture);

current_texture = texture{1};%Choose first texture as initial

Zmin = 0.5;     %min depth
Zmax = 2;       %max depth

%Object init position [m]
Z = Zmax;

%Starting angle - send here a command
cmdstr = sprintf('./checkEnvironment %s %s %d %d left.png right.png %d',current_texture,current_texture,Z,Z,anglenew);
[status, res] = system(cmdstr);   
tmpX = regexp(res,'*','split');        

tic %start time count

while(~status_done) %run until you get the quit signal
  
    if (strcmp(phase,'testing'))
        testing_counter = testing_counter+1; %iteration counter
        %check # of iterations
        if (testing_counter > testing_window)
            phase = 'training'
            testing_counter = 0;
            save(savepath_tests,'model_data')
        end

        % switch object
        if (~mod(testing_counter-1,model.Interval) && testing_counter~=1 && t~=1)
            current_texture = texture{(randi(n_textures,1))};    
            % random position
            Z = Zmin+(Zmax-Zmin)*rand(1,1);% random depth                
            % reset vergence to random value
            anglenew = randi(vergemax,1);% relax the eyes
            cmdstr = sprintf('./checkEnvironment %s %s %d %d left.png right.png %d',current_texture,current_texture,Z,Z,anglenew);
            [status, res] = system(cmdstr);                           
            tmpX = regexp(res,'*','split');                      
        end            
    end
    
    if (strcmp(phase,'training'))
        t = t+1; %iteration counter
        %check # of iterations
        if (t > model.train_time)
        %status_done = 1;
            phase = 'testing'
        end
        
        % switch object
        if ~mod(t-1,model.Interval)
            current_texture = texture{(randi(n_textures,1))};      
            % random position
            Z = Zmin+(Zmax-Zmin)*rand(1,1);% random depth                
            % reset vergence to random value
            anglenew = randi(vergemax,1);% relax the eyes
            cmdstr = sprintf('./checkEnvironment %s %s %d %d left.png right.png %d',current_texture,current_texture,Z,Z,anglenew);
            [status, res] = system(cmdstr);                           
            tmpX = regexp(res,'*','split');              
        end        
    end    
       
    % READ INPUT IMAGES    
    ImageLeft = imread('left.png');
    ImageRight = imread('right.png');                        
    imgrawLeft = ImageLeft;
    imgrawRight = ImageRight;
    h1 = 240;
    w1 = 320;
    h2 = 240;
    w2 = 320;
        
    %LEFT image
        
    img = imgrawLeft;
    %convert to gray scale
    img = .2989*img(:,:,1) +.5870*img(:,:,2) +.1140*img(:,:,3);
    
    
    %SMALL SCALE (16x16)
        
    %Downsample image to 8x8 using Gaussian Pyramid
    img_S = img;
    for i = 1:log2(Dsratio_S)
        img_S = impyramid(img_S,'reduce');
    end
        
    %convert to double
    img_S = double(img_S);
    
    %cut fovea in the center
    [h,w,depth] = size(img_S);
    img_S = img_S(fix(h/2+1-fovea_S/2):fix(h/2+fovea_S/2), fix(w/2+1-fovea_S/2):fix(w/2+fovea_S/2));
        
    %cut patches and store them as col vectors
    patches = im2col(img_S,[patch_size patch_size],'sliding'); %slide window of 1 px
    patches_no = im2col(img_S,[patch_size patch_size],'distinct');    %no overlapping patches (for display)
        
    %take patches at steps of s (8 px)
    cols = [];
    for kc = 1:s_S:nc_S
        C = (kc-1)*nc_S+1:s_S:kc*nc_S;
        cols = [cols C];
    end
        
    patches = patches(:,cols);       %81 patches    

    %pre-processing steps (0 mean, unit norm)
    patches = patches - repmat(mean(patches),[size(patches,1) 1]);  %0 mean
    normp = sqrt(sum(patches.^2));    %patches norm
    %normalize patches to norm 1
    normp(normp==0) = eps;      %regularizer
    patches = patches./repmat(normp,[size(patches,1) 1]); %normalized patches
        
        
    %save patches
    patchesLeft_S = patches;
    patchesLeft_Sno = patches_no;
        
        
    %LARGE SCALE (64x64)
        
    %Downsample image to 8x8 using Gaussian Pyramid
    img_L = img;
    for i = 1:log2(Dsratio_L)
        img_L = impyramid(img_L,'reduce');
    end
        
    %convert to double
    img_L = double(img_L);
        
        
    %cut fovea in the center
    [h,w,depth] = size(img_L);
    img_L = img_L(fix(h/2+1-fovea_L/2):fix(h/2+fovea_L/2), fix(w/2+1-fovea_L/2):fix(w/2+fovea_L/2));
        
    %cut patches and store them as col vectors
    patches = im2col(img_L,[patch_size patch_size],'sliding'); %slide window of 1 px
    patches_no = im2col(img_L,[patch_size patch_size],'distinct');    %no overlapping patches (for display)
        
        
    %take patches at steps of s (8 px)
    cols = [];
    for kc = 1:s_L:nc_L
        C = (kc-1)*nc_L+1:s_L:kc*nc_L;
        cols = [cols C];
    end
        
    patches = patches(:,cols);       %81 patches
        
    patches = patches - repmat(mean(patches),[size(patches,1) 1]);  %0 mean
    normp = sqrt(sum(patches.^2));    %patches norm
        
    %normalize patches to norm 1
    normp(normp==0) = eps;      %regularizer
    patches = patches./repmat(normp,[size(patches,1) 1]); %normalized patches
        
    %store patch at current disparity k
    patchesLeft_L = patches;
    patchesLeft_Lno = patches_no;    
    
    %RIGHT image
    img = imgrawRight;
    %convert to gray scale
    img = .2989*img(:,:,1) +.5870*img(:,:,2) +.1140*img(:,:,3);
           
    %SMALL SCALE (16x16)
        
    %Downsample image to 8x8 using Gaussian Pyramid
    img_S = img;
    for i = 1:log2(Dsratio_S)
        img_S = impyramid(img_S,'reduce');
    end
        
    %convert to double
    img_S = double(img_S);
        
        
    %cut fovea in the center
    [h,w,depth] = size(img_S);
    img_S = img_S(fix(h/2+1-fovea_S/2):fix(h/2+fovea_S/2), fix(w/2+1-fovea_S/2):fix(w/2+fovea_S/2));
        
    %cut patches and store them as col vectors
    patches = im2col(img_S,[patch_size patch_size],'sliding'); %slide window of 1 px
    patches_no = im2col(img_S,[patch_size patch_size],'distinct');    %no overlapping patches (for display)
        
        
    %take patches at steps of s (8 px)
    cols = [];
    for kc = 1:s_S:nc_S
        C = (kc-1)*nc_S+1:s_S:kc*nc_S;
        cols = [cols C];
    end
        
    patches = patches(:,cols);       %81 patches
        
    patches = patches - repmat(mean(patches),[size(patches,1) 1]);  %0 mean
    normp = sqrt(sum(patches.^2));    %patches norm
        
    %normalize patches to norm 1
    normp(normp==0) = eps;      %regularizer
    patches = patches./repmat(normp,[size(patches,1) 1]); %normalized patches
        
    %store patch at current disparity k
    patchesRight_S = patches;
    patchesRight_Sno = patches_no;
        
        
    %LARGE SCALE (64x64)
        
    %Downsample image to 8x8 using Gaussian Pyramid
    img_L = img;
    for i = 1:log2(Dsratio_L)
        img_L = impyramid(img_L,'reduce');
    end
        
    %convert to double
    img_L = double(img_L);
        
        
    %cut fovea in the center
    [h,w,depth] = size(img_L);
    img_L = img_L(fix(h/2+1-fovea_L/2):fix(h/2+fovea_L/2), fix(w/2+1-fovea_L/2):fix(w/2+fovea_L/2));
        
    %cut patches and store them as col vectors
    patches = im2col(img_L,[patch_size patch_size],'sliding'); %slide window of 1 px
    patches_no = im2col(img_L,[patch_size patch_size],'distinct');    %no overlapping patches (for display)
        
        
    %take patches at steps of s (8 px)
    cols = [];
    for kc = 1:s_L:nc_L
        C = (kc-1)*nc_L+1:s_L:kc*nc_L;
        cols = [cols C];
    end
        
    patches = patches(:,cols);       %81 patches
    
    patches = patches - repmat(mean(patches),[size(patches,1) 1]);  %0 mean
    normp = sqrt(sum(patches.^2));    %patches norm

    %normalize patches to norm 1
    normp(normp==0) = eps;      %regularizer
    patches = patches./repmat(normp,[size(patches,1) 1]); %normalized patches
        
    %store patch at current disparity k
    patchesRight_L = patches;
    patchesRight_Lno = patches_no;
            
    %image patches matrix (input to model)
    Current_View = {[patchesLeft_L;patchesRight_L] [patchesLeft_S;patchesRight_S]};
    
    if (strcmp(phase,'training'))
    
        %Compute Vergence Command
        fixdepth = (0.5*baseline)/tand(anglenew/2);
        %save fixaton and object depth    
        model.Z = [model.Z;Z];
        model.fixZ =[model.fixZ;fixdepth];

        %compute current disparity
        angledes = 2*atan(baseline/(2*Z)); %desired vergence [rad]
        disparity = 2*f*tan((angledes-anglenew*pi/180)/2); %current disp [px]      
        model.disp_hist = [model.disp_hist;disparity];

        %compute current vergence error
        anglerr = angledes*180/pi-anglenew; %vergence error [deg]
        %save verge error
        model.vergerr_hist = [model.vergerr_hist;anglerr];

        %generate input feature vector from current images
        [feature,reward,Error,Error_L,Error_S] = model.generateFR(Current_View);
        model.recerr_hist(t,:) = [Error_L;Error_S];
        model.verge_actual = [model.verge_actual;anglenew]; %current vergence angle

        if isempty(model.learned_file)
            %train 2 sparse coding models
            model.scmodel_Large.stepTrain(Current_View{1});
            model.scmodel_Small.stepTrain(Current_View{2});            
            command = model.rlmodel.stepTrain(feature,reward,mod(t-1,model.Interval));    
        end
        
        if (strcmp(phase,'training'))
            sprintf('\n Training Iteration = %.3g \n Command = %.3g  Current Vergence = %.3g Vergence Error = %.3g \n Rec Error = %.3g ',t,command,anglenew,anglerr,Error)
        end
        
        anglenew = max(0.01,command+anglenew); %command is relative angle - constrain to positive vergence
        %safety on control command - RESET TO a new vergence angle
        if anglenew > vergemax
            anglenew = randi(vergemax,1);     %relax the eyes
        end
        cmdstr = sprintf('./checkEnvironment %s %s %d %d left.png right.png %d',current_texture,current_texture,Z,Z,anglenew);    
        % pause(0.1)
        [status, res] = system(cmdstr);    
        tmpX = regexp(res,'*','split');  
        
        %display % completed for training and save model
        if(~mod(t,model.train_time/10) && isempty(model.learned_file) )
            disp([num2str(t/model.train_time*100) '% is finished']);
            save(savepath,'model')

            %save Basis
            model.scmodel_Large.saveBasis;
            model.scmodel_Small.saveBasis;

            %save Weights
            model.rlmodel.saveWeights;  %save policy and value net weights            

            phase = 'testing'
        end
      
    elseif (strcmp(phase,'testing'))
    
        %Compute Vergence Command
        fixdepth = (0.5*baseline)/tand(anglenew/2);
        %save fixaton and object depth    
        model_data.Z = [model_data.Z;Z];
        model_data.fixZ =[model_data.fixZ;fixdepth];

        %compute current disparity
        angledes = 2*atan(baseline/(2*Z)); %desired vergence [rad]
        disparity = 2*f*tan((angledes-anglenew*pi/180)/2); %current disp [px]      
        model_data.disp_hist = [model_data.disp_hist;disparity];

        %compute current vergence error
        anglerr = angledes*180/pi-anglenew; %vergence error [deg]
        %save verge error
        model_data.vergerr_hist = [model_data.vergerr_hist;anglerr];

        %generate input feature vector from current images
        [feature,reward,Error,Error_L,Error_S] = model.generateFR(Current_View);
        model_data.recerr_hist(end+1,:) = [Error_L;Error_S];
        model_data.verge_actual = [model_data.verge_actual;anglenew]; %current vergence angle

        if isempty(model.learned_file)
            command = model.rlmodel.softmaxAct(feature);
        end
        
        if (strcmp(phase,'testing'))
            sprintf('\n Testing Iteration = %.3g \n Command = %.3g  Current Vergence = %.3g Vergence Error = %.3g \n Rec Error = %.3g ',testing_counter,command,anglenew,anglerr,Error)
        end
        
        anglenew = max(0.01,command+anglenew); %command is relative angle - constrain to positive vergence
        %safety on control command - RESET TO a new vergence angle
        if anglenew > vergemax
            anglenew = randi(vergemax,1);     %relax the eyes
        end
        cmdstr = sprintf('./checkEnvironment %s %s %d %d left.png right.png %d',current_texture,current_texture,Z,Z,anglenew);    
        % pause(0.1)
        [status, res] = system(cmdstr);    
        tmpX = regexp(res,'*','split'); 
        
        if ((t >= model.train_time) && (testing_counter+1 > testing_window))
            status_done = 1
        end
        
    end      
end

TotT = toc/60; %total simulation time
sprintf(' Time [min] = %.2f',TotT)

%save model when trained
save(savepath,'model')
save(savepath_tests,'model_data')

end

