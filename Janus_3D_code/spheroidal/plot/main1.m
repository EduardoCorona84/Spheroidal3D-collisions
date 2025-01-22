% python colormap
%
% 05/28/24

clear
addpath('./utils/')

%%% you need to pass your python installation path to matlab, for example:
% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.9.13/bin/python"');
cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/python3"');
% cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/Cellar/python@3.11/3.11.5/bin/python3"');

[xx,yy] = meshgrid(linspace(-1,1,21));
vals = rand(size(xx));

figure(1),clf,
imagesc(vals)
colormap(cmp), caxis([0 1]); colorbar
