mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
%%
new_A_list = {};
new_b_list = {};
for i = 1:5 
    load(fullfile(dirname, ['5x5x5_amphi_data_' num2str(i) '.mat']), ...
        'A_list', 'b_list');
    new_A_list = cat(2, new_A_list,A_list);
    new_b_list = cat(2, new_b_list,b_list);
end
A_list = new_A_list;
b_list = new_b_list;
save(fullfile(dirname, 'amphi_5x5x5.mat'), ...
    'A_list', 'b_list');
%%
figure
mm = min([errStruct_pq.f', errStruct_nic.f', errStruct_stephen.f]) - 1e-10;
semilogy(errStruct_pq.f - mm,'LineWidth',3);
hold on
semilogy(errStruct_nic.f - mm,'LineWidth',3);
semilogy(errStruct_stephen.f' - mm,'--','LineWidth',3);
legend({'ProxQN', '0SR1-Nic', '0SR1-Stephen'})
%%
new_A_list = {};
new_b_list = {};
prefixs = {'bimetallic_2x2x2','bimetallic_3x3x3','amphi_4x4x4', 'amphi_5x5x5'};
for ix = 1:numel(prefixs)
    fname = prefixs{ix};
    load([fname '.mat'], ...
        'A_list', 'b_list');
    new_A_list = cat(2, new_A_list,A_list);
    new_b_list = cat(2, new_b_list,b_list);
end
A_list = new_A_list;
b_list = new_b_list;
save(fullfile(dirname, 'all_data.mat'), ...
    'A_list', 'b_list');