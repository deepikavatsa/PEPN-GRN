function [all_edges] = weighted_edge_prob( prod_con_op_sign, decay_con_op_sign, sus_decay_con_op_sign, sus_prod_con_op_sign, num_of_genes,produce_count,decay_count,sus_decay_count,sus_prod_count)

%%============= GENERATING OCCURENCE MATRIC FOR EACH EVIDENCE ===================

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
%disp('sus_prod_control_edges Done')

%%============ FINDING COMBINED CONTROL EDGES BY COMBINING INFO FROM EACH
%%PROBABILITY MATRIX OF EACH CATEGORY ========================================
%save 'prod_control_edges.mat' prod_control_edges
%save 'decay_control_edges.mat' decay_control_edges
%save 'sus_decay_control_edges.mat' sus_decay_control_edges
%save 'sus_prod_control_edges.mat' sus_prod_control_edges

%% =============== Weighted probability computation ====================

all_edges=[];
for i = 1 : num_of_genes
    output_place = i;
        
    prod_cnt = produce_count(i,1);
    sus_prod_cnt = sus_prod_count(i,1);
    dec_count = decay_count(i,1);
    sus_decay_cnt = sus_decay_count(i,1);
    cnt = prod_cnt + sus_prod_cnt + dec_count + sus_decay_cnt;
        
    all_op1 = prod_edge_occ(:,2);
    all_op2 = sus_prod_edge_occ(:,2);
    all_op3 = decay_edge_occ(:,2);
    all_op4 = sus_decay_edge_occ(:,2);
   
    c1 = prod_edge_occ(all_op1 == output_place,:);
    c11 = c1(: , 1:3);
    c2 = sus_prod_edge_occ(all_op2 == output_place,:);
    c22 = c2(: , 1:3); 
    c3 = decay_edge_occ(all_op3 == output_place,:);
    c33 = c3(: , 1:3);
    c4 = sus_decay_edge_occ(all_op4 == output_place,:);
    c44 = c4(: , 1:3);
    
    edge1 = intersect(c11,c22,'rows');          % common edges of prod and sus prod evidence
    edge2 = intersect(c33,c44,'rows');          % common edges of decay and sus_decay evidence
    edge3 = intersect(edge1,edge2,'rows');      % common edges of all four evidence

    [~,idx1,~] = intersect(c11,edge3,'rows');   % indexes of common edges in prod evidence
    [~,idx2,~] = intersect(c22,edge3,'rows');   % indexes of common edges in sus_prod evidence
    [~,idx3,~] = intersect(c33,edge3,'rows');   % indexes of common edges in decay evidence
    [~,idx4,~] = intersect(c44,edge3,'rows');   % indexes of common edges in sus_decay evidence
    cnt1 = c1(idx1 , 4);                       % cnt of common edges in prod evidence
    cnt2 = c2(idx2 , 4);
    cnt3 = c3(idx3 , 4);                       
    cnt4 = c4(idx4 , 4);
    prob = (cnt1 + cnt2 + cnt3 + cnt4) ./ cnt;
    all_edges = [all_edges ; [edge3 , prob]];
    
    [other_edges1 , ~] = setdiff(c11,edge3 ,'rows');  % edges only present in c11
    [other_edges2 , ~] = setdiff(c22,edge3 ,'rows');  % edges only present in c22
    [other_edges3 , ~] = setdiff(c33,edge3 ,'rows');  % edges only present in c33
    [other_edges4 , ~] = setdiff(c44,edge3 ,'rows');  % edges only present in c44
    other_edges = [other_edges1 ; other_edges2 ; other_edges3 ; other_edges4];
    other_edges(~any(other_edges,2),:) = [];
    other_edges = unique(other_edges,'rows');           % all other edges which are not present in atleast of the evidence edge set
    all_other_edges = zeros(size(other_edges,1) , 4);
    clear cnt1 cnt2 edge prob
    for j = 1 : size(other_edges,1)
       edge = other_edges(j,:);
       [flag1,index1] = ismember(edge,c11,'rows'); 
       [flag2,index2] = ismember(edge,c22,'rows');
       [flag3,index3] = ismember(edge,c33,'rows'); 
       [flag4,index4] = ismember(edge,c44,'rows');

           if flag1 == 1
               cnt1 = c1(index1,4);
           else
               cnt1 = 0;
           end

           if flag2 == 1
               cnt2 = c2(index2,4);
           else
               cnt2 = 0;
           end
           
           if flag3 == 1
               cnt3 = c3(index3,4);
           else
               cnt3 = 0;
           end
           
           if flag4 == 1
               cnt4 = c4(index4,4);
           else
               cnt4 = 0;
           end

           prob = ((cnt1 + cnt2 + cnt3 + cnt4 ) / cnt);        
           all_other_edges(j,:) = [edge , prob];
    end    
    all_edges = [all_edges ; all_other_edges];
    
end



end











