clear
load 'People.mat';

options = pca('options');
options.plots = 'none';
options.preprocessing = {'autoscale'};

cvoptions = crossval('options');
cvoptions.preprocessing = 2;
options.plots = 'none';

ncomp = 4;
model = pca(data, ncomp, options);
modelcv = crossval(data, [], 'pca', {'loo'}, ncomp, cvoptions);

x1 = model.tsqs{1, 1};
y1 = model.ssqresiduals{1, 1};
x2 = modelcv.tsqs{1, 1};
y2 = modelcv.ssqresiduals{1, 1};
figure
hold on
scatter(x1, y1, '.b');
scatter(x2, y2, '.r');
hold off
text(x1, y1, obj_names);
text(x2, y2, obj_names);
line([model.detail.tsqlim{1} model.detail.tsqlim{1}], [min(y1) max(y1)])
line([min(x1) max(x1)], [model.detail.reslim{1} model.detail.reslim{1}])

disp(model.detail.reslim)