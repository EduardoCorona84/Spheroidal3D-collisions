function [ps, gmresTols, out] = zipDenseMats(srcDir)
%% set path
[dirname, ~] = setPaths();
%% set defaults
if ~exist('srcDir', 'var') || isempty(srcDir)
    srcDir = fullfile(dirname, '../data/amphiLCPs.n_2.p_8.cDist_2.3');
end
%% Set all outputs
ps = [];
gmresTols = [];
elems = dir(srcDir);
%%
out = cell(1,500);
for i = 1:numel(elems)
    elem = elems(i);
    if elem.isdir || ~contains(elem.name, '.mat') || contains(elem.name, 'allMats')
        continue
    end
    filePath = fullfile(elem.folder, elem.name);
    prts = split(elem.name, '.mat');
    ix = str2double(prts{1});
    res = load(filePath);
    if isempty(ps) || isempty(gmresTols)
        ps = res.ps;
        gmresTols = res.gmresTols;
    else
        assert(all( ps == res.ps) && all(gmresTols == res.gmresTols));
    end
    out{ix} = res.out;
end
%%
save(fullfile(srcDir, 'allMats.mat'), 'ps', 'gmresTols', 'out');

