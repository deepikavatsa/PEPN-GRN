function [final_edges] = edge_sign_selection(all_edges)

ff = all_edges(:,[1:2]);                            % reg-target pairs from all_edges
[uni,~,id2] = unique(ff,'rows','stable');			% identify unique pairs  in 'uni'
occ = histc(id2,unique(id2));                       % occ of each pair in 'occ'	(Occurrence can be either 1(single regulation sign) or 2(edge with both regulation sign))
idx = find(occ > 1);                                % indexes of edges which occur more than once
duplicate_edges = uni(idx,:);                       % identify edges which occr more than once
single_edges = setdiff(uni,duplicate_edges,'rows');		% rest of the edges occur just once

[~, single_edge_idx] = intersect(ff,single_edges,'rows','stable');	% identify indexes of single_edges in ff
% OR [~,single_edge_idx] = ismember(single_edges,ff,'rows');

dup_edge_idx = setdiff(1 : size(ff,1) , single_edge_idx);	% rest of the indexes in ff represent multiple edges indexes
dup_edge_idx = dup_edge_idx';
edges_with_both_sign = all_edges(dup_edge_idx , :);			% retrieve edges with both signs


aa = edges_with_both_sign(:,[1,2,4]);


% =========== section 2 ===============
[Ua,~,ix] = unique(aa(:,1:2),'rows');                        % Unique Values & Indices
L = accumarray(ix, aa(:,3), [], @(x)max(sum(x)~=2*x));       % Eliminate Duplicates (Returns Logical Vector)
R = accumarray(ix, aa(:,3), [], @max);                       % Maximum Probabilities
C = [Ua(L,:), R(L)] ;
[~,iddx] = intersect(aa,C,'rows','stable');
sel_edges = edges_with_both_sign(iddx,:);
sel_edges2 = all_edges(single_edge_idx,:);
final_edges = [sel_edges ; sel_edges2];

end
