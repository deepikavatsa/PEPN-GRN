function [] = code_v2(num_of_genes,netnum,disc,disclevel,variant)

tic

load (strcat('all_data_',num2str(num_of_genes),'_',num2str(netnum),'_',disc,'_',disclevel,'.mat'))
load (strcat('all_experiment_',num2str(num_of_genes),'_',num2str(netnum),'.mat'))

%% ==========  Compute state pairs  =========================

state_pairs={}; Diff_vec=[];
for i = 1 :size(all_experiment,2)
	
    exp = all_experiment{1,i};
    exp_data = all_data( :,exp(1,1):exp(1,end) );   % extract experimental data
    
    data_pairs = cell(1 , size(exp_data,2)-1);
    for j = 1: (size(exp_data,2)-1)
    data_pairs{j} = [ exp_data(:,j),exp_data(:,j+1) ];
    Diff_vec = [Diff_vec , (exp_data(:,j+1) - exp_data(:,j)) ];     % Difference vector matrix
    end
    state_pairs = [state_pairs,data_pairs];     % cell array of all state pairs
end


%% ====================== PARAMETERS ==========================
num_of_state_pairs = size(state_pairs,2);

allowed_values = [0,1,2];
allowed_changes = [-2,-1,0,1,2];
regulators = [];                                % predefined regulator indexes
if isempty(regulators)  
    num_of_reg = num_of_genes;                  % num_of_reg can be set as sqrt(num_of_genes) for large data sets
else
    num_of_reg = length(regulators);
end


% Pre allocation
prod_con_op_sign = cell(num_of_state_pairs,1);
decay_con_op_sign = cell(num_of_state_pairs,1);
sus_decay_con_op_sign = cell(num_of_state_pairs,1);
sus_prod_con_op_sign = cell(num_of_state_pairs,1);

sus_decay_count = zeros(num_of_genes,1);
sus_prod_count = zeros(num_of_genes,1);

p = cell(1,num_of_state_pairs);
d = cell(1,num_of_state_pairs);
sd = cell(1,num_of_state_pairs);
sp = cell(1,num_of_state_pairs);

%% ================ EVIDENCE COMPUTATION =======================
delete(gcp('nocreate'));
%parpool = 4;

parfor j = 1 : num_of_state_pairs                  % For each state pair
    %disp('State Pair: ');
    %j
    p_str = zeros(num_of_genes,1); 
    d_str = zeros(num_of_genes,1); 
    sd_str = zeros(num_of_genes,1);                % column vector to count sus_decay events in each state pair
    sp_str = zeros(num_of_genes,1);                % column vector to count sus_prod events in each state pair
 
    state_pair = state_pairs{1,j};
    S2 = state_pair(:,2); S1 = state_pair(:,1);
    D = S2 - S1;          			   % Difference vector of this state pair
       
    %%====================== PRODUCTION EVIDENCE ===================================
    prod_op_indexes = find(D > 0);                  % Production evidence in single time difference
            
    if ~isempty(prod_op_indexes)  
        nn = size(prod_op_indexes,1);               	  % nn is num of gene produced
        prod_con = zeros(nn * (num_of_genes-1) , 3);      % pre-allocation. prod_con will contain all edges retrieved from production evidence
                                                  
        for l = 1 : nn
            gene_idx = prod_op_indexes(l,1);
            exp_value_pre = S1(gene_idx);
            exp_value_post = S2(gene_idx);
           
            if exp_value_pre == 1 && exp_value_post == 2 ||   exp_value_pre == 0 && exp_value_post == 2
               p_str(gene_idx,1) =  1;
            [prod_cp,prod_cs] = prod_evidence(S1, gene_idx, regulators);
            k = size(prod_cp,1);
            arr=repmat(gene_idx,k,1);
            edges = [prod_cp,arr,prod_cs];                % concatenation of control, output gene, control sign
                
            lastidx = find(any(prod_con,2),1,'last');     % last filled row in prod_con
            if isempty(lastidx)
                id = 1;
            else
                id = lastidx+1;
            end
            prod_con(id:id+k-1 , :) = edges;       	  % appending newly obtained edges in the prod_con matrix
            
	    end
                
       end
       p{j} = p_str;
       prod_con( ~any(prod_con,2), : ) = [];        	% omit free space(rows with all zeros) in the matrix. 
                                                    
       prod_con_op_sign{j} = prod_con;              
    end
	%%============= DECAY EVIDENCE =============================
    decay_op_indexes = find(D < 0);                 	% Decay evidence in single time difference
            
    if ~isempty(decay_op_indexes)  
        nn = size(decay_op_indexes,1);              	% nn is num of genes decayed
        decay_con = zeros(nn * (num_of_genes-1) , 3);     % pre-allocation. decay_con will contain all edges retrieved from decay evidence
        for l = 1 : nn
            gene_idx = decay_op_indexes(l,1);
            exp_value_pre = S1(gene_idx);
            exp_value_post = S2(gene_idx);
           
            if exp_value_pre == 1 && exp_value_post == 0 || exp_value_pre == 2 && exp_value_post == 0
                d_str(gene_idx,1) =  1;
            [decay_cp,decay_cs] = decay_evidence(S1, gene_idx, regulators);
            
            k = size(decay_cp,1);
            arr=repmat(gene_idx,k,1);
            edges = [decay_cp,arr,decay_cs];         	% concatenation of control, output gene, control sign
            
            lastidx = find(any(decay_con,2),1,'last');
            if isempty(lastidx)
                id = 1;
            else
                id = lastidx+1;
            end
            decay_con(id:id+k-1 , :) = edges;       	% appending newly obtained edges in the decay_con matrix

	    end
        end
	d{j} = d_str;
        decay_con( ~any(decay_con,2), : ) = [];     	% omit free space(rows with all zeros) in the matrix
                                                    
        decay_con_op_sign{j} = decay_con;           
            
    end
    
    %%================= ADDITIONAL EVIDENCE ========================
    gene_idxes = 1 : num_of_genes;
        gene_idxes([prod_op_indexes',decay_op_indexes']) = [];  % genes left after filtering out produced or decayed genes in the difference vector
           
        nn = size(gene_idxes,2);                            	% nn is genes that can get sus_decay or sus_prod in state pair
        sus_decay_con = zeros(nn * (num_of_genes-1) , 3);   	% pre-allocation
        sus_prod_con = zeros(nn * (num_of_genes-1) , 3);    	% pre-allocation
        for l = 1 : nn
            gene_idx = gene_idxes(1,l);
            exp_value_pre = S1(gene_idx);
            exp_value_post = S2(gene_idx);
           
            if exp_value_pre == 0 && exp_value_post == 0    	% sustained decay evidence
                sd_str(gene_idx,1) =  1;
                
                [sd_cp, sd_cs] = sus_decay_evidence(S1, gene_idx, regulators);
                       
                k = size(sd_cp,1);
                arr=repmat(gene_idx,k,1);
                edges = [sd_cp,arr,sd_cs];
                
                lastidx = find(any(sus_decay_con,2),1,'last');
                if isempty(lastidx)
                    id = 1;
                else
                    id = lastidx+1;
                end
                sus_decay_con(id:id+k-1 , :) = edges;
            elseif exp_value_pre == 2 && exp_value_post == 2 	% sustained production evidence
                sp_str(gene_idx,1) = 1;
                
                [sp_cp, sp_cs] = sus_prod_evidence(S1, gene_idx, regulators);
                
                k = size(sp_cp,1);
                arr=repmat(gene_idx,k,1);
                edges = [sp_cp,arr,sp_cs];
                
                lastidx = find(any(sus_prod_con,2),1,'last');
                if isempty(lastidx)
                    id = 1;
                else
                    id = lastidx+1;
                end
                sus_prod_con(id:id+k-1 , :) = edges;
            end
        end
        sd{j} = sd_str;
        sp{j} = sp_str;
        
        sus_decay_con( ~any(sus_decay_con,2), : ) = [];     % omit free space(rows with all zeros) in the matrix
        sus_prod_con( ~any(sus_prod_con,2), : ) = [];       % omit free space(rows with all zeros) in the matrix
        
        sus_decay_con_op_sign{j} = sus_decay_con;       	% sus decay con edges for this state pair
        sus_prod_con_op_sign{j} = sus_prod_con;         	% sus prod con edges for this state pair
        
end

delete(gcp('nocreate'));

d = cell2mat(d);
decay_count = sum(d , 2);

p = cell2mat(p);
produce_count = sum(p , 2);

sd = cell2mat(sd);
sus_decay_count = sum(sd , 2);

sp = cell2mat(sp);
sus_prod_count = sum(sp , 2);

prod_con_op_sign = cell2mat(prod_con_op_sign);
decay_con_op_sign = cell2mat(decay_con_op_sign);
sus_prod_con_op_sign = cell2mat(sus_prod_con_op_sign);
sus_decay_con_op_sign = cell2mat(sus_decay_con_op_sign);

for i = 1:num_of_genes
   occ_count(i,1) = size(find(all_data(i,:) > 0),2);
   decay_edges(i,1) = decay_count(i,1)/occ_count(i,1);
end

%% Probability computation of each control edge. Also, decay probability is also computed for each gene.

[all_edges] = weighted_edge_prob(prod_con_op_sign,decay_con_op_sign,sus_decay_con_op_sign,sus_prod_con_op_sign,num_of_genes,produce_count,decay_count,sus_decay_count,sus_prod_count);


[final_edges] = edge_sign_selection(all_edges);

% ====================== Edge file ==============================
score = final_edges(:,[1,2,4]);
[~,index] = sort(score(:,3),'descend');

load (strcat('genenames_',num2str(num_of_genes),'gene.mat'))

filename=strcat('predicted_edges_',num2str(num_of_genes),'_',num2str(netnum),'_',disc,'.txt');
fileID = fopen(filename,'w');
for i = 1:size(index,1)
   idx = index(i,1);
   reg = score(idx,1);
   reg_name = genenames{reg,1};
   target = score(idx,2);
   target_name = genenames{target,1};
   prob = score(idx,3);
   fprintf(fileID,'%s\t%s\t%s\n',reg_name,target_name,num2str(prob));
end

fclose(fileID);

% ====================== Probability Matrix =====================

num = num_of_genes*num_of_genes;

prob_matrix = zeros(num_of_genes , num_of_genes);
prob_matrix(logical(eye(size(prob_matrix)))) = NaN;

prob_sign = NaN(num_of_genes , num_of_genes);

for i = 1 : size(final_edges,1)
   reg_idx = final_edges(i,1);
   tar_idx = final_edges(i,2);
   sign = final_edges(i,3);
   prob = final_edges(i,4);
   prob_matrix(tar_idx , reg_idx) = prob;
   prob_sign(tar_idx , reg_idx) = sign;
end

prob_matrix = prob_matrix';
prob_sign = prob_sign';
probability = reshape(prob_matrix,num,1);
prob_sign = reshape(prob_sign,num,1);

pr = strcat('prob_', num2str(num_of_genes), '_' , num2str(netnum),'.mat');
pr_sign = strcat('prob_sign_', num2str(num_of_genes), '_' , num2str(netnum),'.mat');
save(pr, 'probability');
save(pr_sign, 'prob_sign');

%% ============== FINDING RESTRICTED REGULATOR EDGES (for large data sets) ====================
%% This code selects top nreg regulators for each target gene. If there exists regulators with same probability and
%% nreg somehow selects only few of them, it would give inconsistent results. Therefore, we compute the 
%% prob of last edge in the selected edges (after applying nreg) and find all edges with atleast
%% this prob. Thus, in this case, number of edges selected will be more than nreg.
%% 

nreg = 3;   % number of regulators allowed for each target gene

tar_genes = final_edges(:,2);
ncol = size(final_edges , 2);
res_final_edges = zeros(num_of_genes*num_of_genes , ncol);
j=1;
for i = 1 : num_of_genes
   tar = i;
   
   tar_edges =  final_edges(tar_genes == tar , :);
   if ~isempty(tar_edges)
        tar_edges = sortrows(tar_edges , 4);
        tar_edges = flipud(tar_edges);
         nrows = size(tar_edges , 1);
         if nrows < nreg
             res_final_edges(j : j+nrows-1 , :) = tar_edges;
             j = j+nrows;
         else
             sel_tar_edges = tar_edges(1:nreg , :);
             prob = sel_tar_edges(end,4);
             a = tar_edges(tar_edges(:,4) >= prob , :);
             arow = size(a,1);
             res_final_edges(j:j+arow-1 , :) = a;
             j = j+arow;
         end
   end
end

res_final_edges( ~any(res_final_edges,2), : ) = [];

% ====================== Edge file ==============================
score = res_final_edges(:,[1,2,4]);
[~,index] = sort(score(:,3),'descend');

filename=strcat('res_reg_edges.txt');
fileID = fopen(filename,'w');
for i = 1:size(index,1)
   idx = index(i,1);
   reg = score(idx,1);
   reg_name = genenames{reg,1};
   target = score(idx,2);
   target_name = genenames{target,1};
   prob = score(idx,3);
   fprintf(fileID,'%s\t%s\t%s\n',reg_name,target_name,num2str(prob));
end

fclose(fileID);
 

toc


end






