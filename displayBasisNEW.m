%%% Display the basis according to the binocularity
%TODO test this updated version with model that kept track of their basis
close all
r = 16; c = 18; %how to arrange the basis in rows and cols

% len = size(model.scModel_Small.basisHist,3);  %# of trials saved
numScales = length(model.scModel);
len = size(model.scModel{1}.basisHist, 3);  %# of trials saved
basisTrack = cell(numScales, len);           %variable to store all the saved basis

% copy into the new var
% for j = 1:len
%     basisTrack{1,j} = model.scModel_Large.basisHist(:,:,j);
%     basisTrack{2,j} = model.scModel_Small.basisHist(:,:,j);
% end
for scale = 1:numScales
    if len == 1
        basisTrack{scale, 1} = model.scModel{scale}.basis;
    else
        for j = 1:len
            basisTrack{scale, j} = model.scModel{scale}.basisHist(:, :, j);
        end
    end
end


% k = 11;

% basisTrack = this.drecord.basisTrack(1:len);
% basisTrack{1} = model.scModel_Large.basisHist(:,:,1);
% basisTrack{1} = model.scModel_Small.basisHist(:,:,k);
% checkPoint = 1;

h = figure(1);
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);

% loop over scales
for s = 1:numScales
    % sort basis according to left energy norm
    endBasis = basisTrack{s,end}(1:end/2,:);
    leftEnergy = abs(sum(endBasis.^2)-0.5);
    [~,I] = sort(leftEnergy);

    subplot(1,2,s);
    [di,num] = size(basisTrack{s,1});

    fun1 = @(blc_struct) padarray(padarray(reshape(permute(padarray(reshape(blc_struct.data, sqrt(di / 2), ...
             sqrt(di / 2), 2), [1, 1], 'pre'), [1, 3, 2]), (sqrt(di / 2) + 1) * 2, sqrt(di / 2) + 1), [1, 1], ...
             'post') - 1, [1 1], 'pre') + 1;

    for j = 1:len
        A = basisTrack{s,j}(:,I);
%         B = reshape(A,di*r,c);
        B = reshape(A,di * r,num / r); %hotfix!
        B = B/max(max(abs(B))) + 0.5;
        C = padarray(padarray(blockproc(B,[di,1],fun1)-1,[1 1],'post')+1,[2,2]);
        imshow(C);
        % title(num2str(model.trainTime*0.1*(j-1)));
        if s == 1
            title('Coarse Scale')
        elseif s == 2
            title('Fine Scale')
        end
        drawnow; pause(.01);
    end
end

% Display Histogram of Norm or Energy of Left (non blurred) Bases
for s = 1:numScales
    endBasis = basisTrack{s,end}(1:end/2,:);
    leftEnergy = abs(sum(endBasis.^2));
    figure(2);
    subplot(1,2,s); hold on
    hist(leftEnergy,0:0.1:1); xlabel('Basis norm'); ylabel('# of bases');
    axis([-0.2 1 0 250])
    if s == 1
        title('Coarse Scale')
    else
        title('Fine Scale')
    end
end

% Display Histogram of Binocularity
for s = 1:numScales
    Lefteye = basisTrack{s,end}(1:end/2,:);
    Righteye = basisTrack{s,end}(end/2+1:end,:);
    % binocular index
    b = (sqrt(sum(Lefteye .^ 2)) - sqrt(sum(Righteye .^ 2))) ./ (sqrt(sum(Lefteye .^ 2)) + sqrt(sum(Righteye .^ 2)));
    bins = [-0.85 -0.5 -0.15 0 0.15 0.5 0.85];

    figure(3);
    subplot(1,2,s); hold on
    hist(b,bins)
    n_elements = histc(b,bins); xlabel('Binocularity'); ylabel('# of bases');
    % h = findobj(gca,'Type','patch');
    % set(h,'FaceColor','b','EdgeColor','w')
    if s == 1
        title('Coarse Scale')
    else
        title('Fine Scale')
    end
end

% hold on;
% xlabel('iterations');
% ylabel('verg err');
% temp = [];
% for i = 1:10:length(model.vergerr_hist)
%    temp(end+1) =  model.vergerr_hist(i+8,1);
% end
% plot(temp);
