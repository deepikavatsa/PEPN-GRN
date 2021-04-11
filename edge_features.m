function [] = edge_features(num_of_genes , netnum)

load genepairs.mat
load p_pr.mat
load d_pr.mat
load sp_pr.mat
load sd_pr.mat

sign = genepairs(:,3);
act_pairs = genepairs(sign==1,:);
inh_pairs = genepairs(sign==0,:);

[~,i1] = ismember(act_pairs,genepairs,'rows');
[~,i2] = ismember(inh_pairs,genepairs,'rows');

act_p_pr = p_pr(i1,:);
inh_p_pr = p_pr(i2,:);

act_d_pr = d_pr(i1,:);
inh_d_pr = d_pr(i2,:);

act_sp_pr = sp_pr(i1,:);
inh_sp_pr = sp_pr(i2,:);

act_sd_pr = sd_pr(i1,:);
inh_sd_pr = sd_pr(i2,:);

all_act_edges = [act_pairs, act_p_pr, act_d_pr, act_sp_pr, act_sd_pr];
all_inh_edges = [inh_pairs, inh_p_pr, inh_d_pr, inh_sp_pr, inh_sd_pr];

%-------------------------------------------------------
groundtruth = readtext('groundtruth_edges_signed.tsv','\t');

%% strip off the first character (G), convert to double
gd_pairs = groundtruth(:,1:2);
B = char(gd_pairs{:,1});
D = B(: , 2:end);
gd_regid = str2num(D);

B = char(gd_pairs{:,2});
D = B(: , 2:end);
gd_tarid = str2num(D);

gd_pairs_id = [gd_regid , gd_tarid];

%% convert sign from - to 0 and + to 1 and then to double
gd_sign = groundtruth(:,3);

logicalidx = cellfun(@(x)isequal(x,'-'),gd_sign);
[row,~] = find(logicalidx);

logicalidx = cellfun(@(x)isequal(x,'+'),gd_sign);
[row1,~] = find(logicalidx);

gd_sign = cell2mat(gd_sign);
gd_sign(row,1) = 0;
gd_sign(row1,1) = 1;
gd_sign = double(gd_sign);

gdtruth_edges = [gd_pairs_id , gd_sign];
gdtruth_act_edges = gdtruth_edges(gdtruth_edges(:,3) == 1 ,:);
gdtruth_inh_edges = gdtruth_edges(gdtruth_edges(:,3) == 0 ,:);

%-----------------------------------------------------------
%% pos_act_edges are activatory edges actually present in gdtruth network 
%% and neg_act_edges are activatory edges not present in gdtruth network

[~, idx] = ismember(gdtruth_act_edges , all_act_edges(:,1:3),'rows');
rowidx = 1:size(all_act_edges,1);
rem_idx = setdiff(rowidx , idx);
pos_act_edges = all_act_edges(idx,:);
neg_act_edges = all_act_edges(rem_idx,:);

[~, idx] = ismember(gdtruth_inh_edges , all_inh_edges(:,1:3),'rows');
rowidx = 1:size(all_inh_edges,1);
rem_idx = setdiff(rowidx , idx);
pos_inh_edges = all_inh_edges(idx,:);
neg_inh_edges = all_inh_edges(rem_idx,:);

file1 = strcat('pos_act_edges_', num2str(num_of_genes),'_', num2str(netnum),'.mat');
file2 = strcat('neg_act_edges_', num2str(num_of_genes),'_', num2str(netnum),'.mat');
file3 = strcat('pos_inh_edges_', num2str(num_of_genes),'_', num2str(netnum),'.mat');
file4 = strcat('neg_inh_edges_', num2str(num_of_genes),'_', num2str(netnum),'.mat');

save(file1, 'pos_act_edges');     
save(file2, 'neg_act_edges');
save(file3, 'pos_inh_edges');
save(file4, 'neg_inh_edges');

end
