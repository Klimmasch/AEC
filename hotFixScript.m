% parent = '/home/klimmasch/projects/results/Discount Factor vs Interval/'
% files = dir('/home/klimmasch/projects/results/Discount Factor vs Interval/*00*');
folder = 'Discount Factor vs Interval'
% folder = 'Regularizer_vs_ActorLR'

parent = strcat('/home/aecgroup/aecdata/Results/', folder);
% parent = strcat('/home/klimmasch/projects/results/', folder);
files = dir(sprintf('%s/*00*', parent));


for f = 1:length(files)
    savePath = sprintf('%s/%s', parent, files(f).name)
    model = load(strcat(savePath, '/model.mat'));

    model = model.model;
    model.savePath = savePath;
    if length(model.rlModel.actorLearningRange) == 1
        val = model.rlModel.actorLearningRange(1);
        model.rlModel.actorLearningRange = [val, val];
    end
%     savePathNew = sprintf('%s/16-10-19_2000000iter_1_cluster_regul_%1.0e_actorLR_[%1.2f-%1.2f]', parent, model.rlModel.CActor.regularizer, model.rlModel.actorLearningRange(1), model.rlModel.actorLearningRange(2))
%     model.savePath = savePathNew;
    save(strcat(model.savePath, '/model'), 'model');
end