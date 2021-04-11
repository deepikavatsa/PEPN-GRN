function [all_edges] = unwt_edge_prob( prod_con_op_sign, decay_con_op_sign, sus_decay_con_op_sign, sus_prod_con_op_sign, num_of_genes,produce_count,decay_count,sus_decay_count,sus_prod_count)


%%============= GENERATING OCCURENCE MATRIX FOR EACH EVIDENCE ===================

[uniqueRow,~,uid] = unique(prod_con_op_sign, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
occurrences = histc(uid, unique(uid));
prod_edge_occ = [uniqueRow,occurrences];

%disp('prod_edge_occ Done')

%-------------------------------------------------------------------
[uniqueRow,~,uid] = unique(decay_con_op_sign, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
occurrences = histc(uid, unique(uid));
decay_edge_occ = [uniqueRow,occurrences];

%disp('decay_edge_occ Done')
%-------------------------------------------------------------------
[uniqueRow,~,uid] = unique(sus_decay_con_op_sign, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
occurrences = histc(uid, unique(uid));
sus_decay_edge_occ = [uniqueRow,occurrences];

%disp('sus_decay_edge_occ Done')
%-------------------------------------------------------------------
[uniqueRow,~,uid] = unique(sus_prod_con_op_sign, 'rows', 'stable');           % Unique Values By Row, Retaining Original Order
occurrences = histc(uid, unique(uid));
sus_prod_edge_occ = [uniqueRow,occurrences];

%disp('sus_prod_edge_occ Done')

%%================== GENERATE CONTROL EDGES FOR EACH EVIDENCE =============
[m , n] = size(prod_edge_occ);
prod_control_edges = zeros( m, n );
all_op = prod_edge_occ(:,2);
for i = 1 : num_of_genes
    cnt = produce_count(i,1);
    output_place = i;
    c = prod_edge_occ(all_op == output_place,:);
    k = size(c,1);
    edge_occ_vector = c(:,4);
    edge_prob_vector = edge_occ_vector/cnt;
    edges = [c(:,1:3),edge_prob_vector];
   
    lastidx = find(any(prod_control_edges,2),1,'last');   % last filled row in prod_con
    if isempty(lastidx)
        id = 1;
    else
        id = lastidx+1;
    end
    prod_control_edges(id:id+k-1 , :) = edges;  
end

prod_control_edges(~any(prod_control_edges,2),:) = [];
%disp('prod_control_edges Done')
%--------------------------------------------------------------------
[m , n] = size(decay_edge_occ);
decay_control_edges = zeros( m, n );
all_op = decay_edge_occ(:,2);
for i = 1:num_of_genes
    cnt = decay_count(i,1);
    output_place = i;
    c = decay_edge_occ(all_op == output_place,:);
    k = size(c,1);
    edge_occ_vector = c(:,4);
    edge_prob_vector = edge_occ_vector/cnt;
    edges = [c(:,1:3),edge_prob_vector];
    
    lastidx = find(any(decay_control_edges,2),1,'last');   % last filled row 
    if isempty(lastidx)
        id = 1;
    else
        id = lastidx+1;
    end
    decay_control_edges(id:id+k-1 , :) = edges;  
end

decay_control_edges(~any(decay_control_edges,2),:) = [];
%disp('decay_control_edges Done')
%-----------------------------------------------------------------------
[m , n] = size(sus_decay_edge_occ);
sus_decay_control_edges = zeros( m, n );

all_op = sus_decay_edge_occ(:,2);
for i = 1:num_of_genes
    cnt = sus_decay_count(i,1);
    output_place = i;
    c = sus_decay_edge_occ(all_op == output_place,:);
    k = size(c,1);
    edge_occ_vector = c(:,4);
    edge_prob_vector = edge_occ_vector/cnt;
    edges = [c(:,1:3),edge_prob_vector];
    
    lastidx = find(any(sus_decay_control_edges,2),1,'last');   % last filled row 
    if isempty(lastidx)
        id = 1;
    else
        id = lastidx+1;
    end
    sus_decay_control_edges(id:id+k-1 , :) = edges;  
end

sus_decay_control_edges(~any(sus_decay_control_edges,2),:) = [];
%disp('sus_decay_control_edges Done')
%-----------------------------------------------------------------------
[m , n] = size(sus_prod_edge_occ);
sus_prod_control_edges = zeros( m, n );

all_op = sus_prod_edge_occ(:,2);
for i = 1 : num_of_genes
    cnt = sus_prod_count(i,1);
    output_place = i;
    c = sus_prod_edge_occ(all_op == output_place,:);
    k = size(c,1);
    edge_occ_vector = c(:,4);
    edge_prob_vector = edge_occ_vector/cnt;
    edges = [c(:,1:3),edge_prob_vector];
    
    lastidx = find(any(sus_prod_control_edges,2),1,'last');   % last filled row 
    if isempty(lastidx)
        id = 1;
    else
        id = lastidx+1;
    end
    sus_prod_control_edges(id:id+k-1 , :) = edges;  
end

sus_prod_control_edges(~any(sus_prod_control_edges,2),:) = [];
%disp('sus_prod_control_edges Done')

%%============ FINDING COMBINED CONTROL EDGES BY COMBINING INFO FROM EACH
%%PROBABILITY MATRIX OF EACH CATEGORY ========================================
% save 'prod_control_edges.mat' prod_control_edges
% save 'decay_control_edges.mat' decay_control_edges
% save 'sus_decay_control_edges.mat' sus_decay_control_edges
% save 'sus_prod_control_edges.mat' sus_prod_control_edges

c1 = prod_control_edges;
c11 = c1(:,1:3);
c2 = decay_control_edges;
c22 = c2(:,1:3);
c3 = sus_prod_control_edges;
c33 = c3(:,1:3);
c4 = sus_decay_control_edges;
c44 = c4(:,1:3);
%------------------------------------
    
edge1 = intersect(c11,c22,'rows');          % common edges of prod and decay evidence
edge2 = intersect(c33,c44,'rows');          % common edges of sus_prod and sus_decay evidence
edge3 = intersect(edge1,edge2,'rows');      % common edges of all four evidence

%-------------------------------------------------
[other_edges1 , ~] = setdiff(c11,edge3,'rows');  % edges only present in c11
[other_edges2 , ~] = setdiff(c22,edge3,'rows');  % edges only present in c22
[other_edges3 , ~] = setdiff(c33,edge3,'rows');  % edges only present in c33
[other_edges4 , ~] = setdiff(c44,edge3,'rows');  % edges only present in c44

other_edges = [other_edges1 ; other_edges2 ; other_edges3 ; other_edges4 ];
other_edges = unique(other_edges,'rows');           % all other edges which are not present in atleast of the evidence edge set

other_edges(~any(other_edges,2),:) = [];

rest_edges1 = setdiff(other_edges , other_edges1, 'rows');  % rest_edges1 are edges present in other evidence but not in c1
m = size(rest_edges1,1);
pr = zeros(m,1);
rest_edges_c1 = [rest_edges1 , pr];
c1 = [c1 ; rest_edges_c1];

rest_edges2 = setdiff(other_edges , other_edges2, 'rows');
m = size(rest_edges2,1);
pr = zeros(m,1);
rest_edges_c2 = [rest_edges2 , pr];
c2 = [c2 ; rest_edges_c2];

rest_edges3 = setdiff(other_edges , other_edges3, 'rows');
m = size(rest_edges3,1);
pr = zeros(m,1);
rest_edges_c3 = [rest_edges3 , pr];
c3 = [c3 ; rest_edges_c3];

rest_edges4 = setdiff(other_edges , other_edges4, 'rows');
m = size(rest_edges4,1);
pr = zeros(m,1);
rest_edges_c4 = [rest_edges4 , pr];
c4 = [c4 ; rest_edges_c4];

% After the above step, each matrix c1, c2, c3 and c4 contains same num of
% rows, so extract edges from the matrix

edges = c1(:,1:3);

[~,idx1] = ismember(edges,c1(:,1:3),'rows');   % indexes of common edges in prod evidence
[~,idx2] = ismember(edges,c2(:,1:3),'rows');   % indexes of common edges in decay evidence
[~,idx3] = ismember(edges,c3(:,1:3),'rows');   % indexes of common edges in sus_prod evidence
[~,idx4] = ismember(edges,c4(:,1:3),'rows');   % indexes of common edges in sus_decay evidence
    
prob1 = c1(idx1 , 4);                       % prob of common edges in prod evidence
prob2 = c2(idx2 , 4);
prob3 = c3(idx3 , 4);
prob4 = c4(idx4 , 4);
    
prob = (prob1 + prob2 + prob3 + prob4) ./ 4;    % final prob column vector

all_edges = [edges , prob];
%disp('all edges done')

%-------------------------------------------------------------------
edges_evidence_prob = [edges , prob1, prob2, prob3, prob4];
save 'edges_evidence_prob.mat' edges_evidence_prob

genepairs = edges;
p_pr = prob1;
d_pr = prob2;
sp_pr = prob3;
sd_pr = prob4;

save 'genepairs.mat' genepairs
save 'p_pr.mat' p_pr
save 'd_pr.mat' d_pr
save 'sp_pr.mat' sp_pr
save 'sd_pr.mat' sd_pr
    

end











