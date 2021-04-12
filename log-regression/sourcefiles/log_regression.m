function [] = log_regression(num_of_genes,totalnets)

pos_act_edges_net = cell(totalnets,1);
neg_act_edges_net = cell(totalnets,1);

pos_inh_edges_net = cell(totalnets,1);
neg_inh_edges_net = cell(totalnets,1);


for i = 1 : totalnets
	file = strcat('pos_act_edges_',num2str(num_of_genes),'_',num2str(i),'.mat');
	load (file);
	pos_act_edges_net{i} =  pos_act_edges;
end

for i = 1 : totalnets
	file = strcat('neg_act_edges_',num2str(num_of_genes),'_',num2str(i),'.mat');
	load (file);
	neg_act_edges_net{i} =  neg_act_edges;
end

for i = 1 : totalnets
	file = strcat('pos_inh_edges_',num2str(num_of_genes),'_',num2str(i),'.mat');
	load (file);
	pos_inh_edges_net{i} =  pos_inh_edges;
end

for i = 1 : totalnets
	file = strcat('neg_inh_edges_',num2str(num_of_genes),'_',num2str(i),'.mat');
	load (file);
	neg_inh_edges_net{i} =  neg_inh_edges;
end


numfolds = totalnets;                       % numfolds is same as number of networks we have
act_weights = cell(numfolds,1); 
inh_weights = cell(numfolds,1); 
act_wt_influence = cell(numfolds,1);
inh_wt_influence = cell(numfolds,1); 

train_features_act = cell(numfolds,1);      train_features_inh = cell(numfolds,1);    
test_features_act = cell(numfolds,1);       test_features_inh = cell(numfolds,1); 
train_labels_act = cell(numfolds,1);        train_labels_inh = cell(numfolds,1); 
test_labels_act = cell(numfolds,1);         test_labels_inh = cell(numfolds,1); 
train_scores_act = cell(numfolds,1);       train_scores_inh = cell(numfolds,1); 
test_scores_act = cell(numfolds,1);        test_scores_inh = cell(numfolds,1); 
test_edges_act = cell(totalnets,1);
test_edges_inh = cell(totalnets,1);

for i = 1 : numfolds
    total_nets = 1 : totalnets;
    test_net = i;
    train_net = setdiff(total_nets, test_net);
    
    %--------------  ACTIVATION EDGES  -----------------------
    % EXTRACT TRAIN FEATURES
    pos_features =[]; neg_features =[]; 
    for j = 1 : size(train_net,2)
        netnum = train_net(j);
        p = pos_act_edges_net{netnum};
        pos_features = [pos_features ; p(:,4:end)];
        n = neg_act_edges_net{netnum};
        neg_features = [neg_features ; n(:,4:end)];
    end
        
    pos_labels = ones(size(pos_features,1) , 1);
    neg_labels = zeros(size(neg_features,1) , 1);

    train_features_for_act_edges = [pos_features ; neg_features];
    train_labels_for_act_edges = [pos_labels ; neg_labels];

    % EXTRACT TEST FEATURES
    p = pos_act_edges_net{test_net};
    n = neg_act_edges_net{test_net};
       
    test_pos_labels = ones(size(p,1) , 1);
    test_neg_labels = zeros(size(n,1) , 1);
    
    test_features_for_act_edges = [p(:,4:end) ; n(:,4:end)]; 
    test_labels_for_act_edges = [test_pos_labels ; test_neg_labels];
    
    act_edges = [p(:,1:3) ; n(:,1:3)];
        
    %-------------  INHIBITORY EDGES  -------------------
    % EXTRACT TRAIN FEATURES
    pos_features =[]; neg_features =[]; 
    for j = 1 : size(train_net,2)
        netnum = train_net(j);
        p = pos_inh_edges_net{netnum};
        pos_features = [pos_features ; p(:,4:end)];
        n = neg_inh_edges_net{netnum};
        neg_features = [neg_features ; n(:,4:end)];
    end
        
    pos_labels = ones(size(pos_features,1) , 1);
    neg_labels = zeros(size(neg_features,1) , 1);

    train_features_for_inh_edges = [pos_features ; neg_features];
    train_labels_for_inh_edges = [pos_labels ; neg_labels];

    % EXTRACT TEST FEATURES
    p = pos_inh_edges_net{test_net};
    n = neg_inh_edges_net{test_net};
    
    test_pos_labels = ones(size(p,1) , 1);
    test_neg_labels = zeros(size(n,1) , 1);
    
    test_features_for_inh_edges = [p(:,4:end) ; n(:,4:end)]; 
    test_labels_for_inh_edges = [test_pos_labels ; test_neg_labels];
    
    inh_edges = [p(:,1:3) ; n(:,1:3)];
    
    %-----------------------------------------------------        
    train_features_act{i} = train_features_for_act_edges;
    train_features_inh{i} = train_features_for_inh_edges;
    
    test_features_act{i} = test_features_for_act_edges;
    test_features_inh{i} = test_features_for_inh_edges;
    
    train_labels_act{i} = train_labels_for_act_edges;
    train_labels_inh{i} = train_labels_for_inh_edges;
    
    test_labels_act{i} = test_labels_for_act_edges;
    test_labels_inh{i} = test_labels_for_inh_edges;
    
    test_edges_act{i} = act_edges;
    test_edges_inh{i} = inh_edges;
    
end

% Weight computation using logistic regression
for i = 1 : numfolds
    
    % =========== ACTIVATION EDGES ================
    tr_features = train_features_act{i,1};
    tr_labels = train_labels_act{i,1};
        
    % make balanced data set using SMOTE upsampling
    [SMOTE_features, SMOTE_labels] = smote(tr_features , tr_labels);
    
    % LEARN WEIGHTS OF LOGISTIC REGRESSION MODEL USING TRAIN FEATURES AND LABELS
    [b,dev,stats] = glmfit( SMOTE_features , SMOTE_labels , 'binomial','logit' );
    act_weights{i} = b;
    act_wt_influence{i} = std(SMOTE_features) .* b(2:end,1)';   % b(2:end,1)' is a row vector of weights excluding the intercept
           
    % FIT MODEL ON TRAIN FEATURES
    train_act_score = glmval(b,tr_features,'logit');
    
    ts_features = test_features_act{i,1};
    % FIT MODEL ON TEST FEATURES
    test_act_score = glmval( b,ts_features,'logit' );
    
    % =========== INHIBITION EDGES ================
    tr_features = train_features_inh{i,1};
    tr_labels = train_labels_inh{i,1};
    
    % make balanced data set using SMOTE upsampling
    [SMOTE_features, SMOTE_labels] = smote(tr_features , tr_labels);
    
    % LEARN WEIGHTS OF LOGISTIC REGRESSION MODEL USING TRAIN FEATURES AND LABELS
    [b,dev,stats] = glmfit( SMOTE_features , SMOTE_labels , 'binomial','logit' );
    inh_weights{i} = b;
    inh_wt_influence{i} = std(SMOTE_features) .* b(2:end,1)';   % b(2:end,1)' is a row vector of weights excluding the intercept
            
    % FIT MODEL ON TRAIN FEATURES
    train_inh_score = glmval( b,tr_features,'logit' );
    
    ts_features = test_features_inh{i,1};
    % FIT MODEL ON TEST FEATURES
    test_inh_score = glmval( b,ts_features,'logit' );
    
    %==================================================
    train_scores_act{i} = train_act_score;
    train_scores_inh{i} = train_inh_score;
    
    test_scores_act{i} = test_act_score;
    test_scores_inh{i} = test_inh_score;
end

%% ==========  WRITE WEIGHTS TO FILE  ===========
fid = fopen('performance-evaluation.txt','wt');
fprintf(fid,'\n %s', 'Learned Weights (w0, w1, w2, w3, w4)');
    fprintf(fid,'\n\n%s\n','============ Activation edges ===========');
    
    for i = 1 : numfolds
	fprintf(fid,'%s%d\t','Fold ',i);
        fprintf(fid,'%s \t %f \t %f \t %f \t %f \t %f \n', 'Weights: ', act_weights{i}' );
    end
    
    fprintf(fid,'\n\n%s\n','============ Inhibition edges ===========');
    for i = 1 : numfolds
	fprintf(fid,'%s%d\t','Fold ',i);
        fprintf(fid,'%s \t %f \t %f \t %f \t %f \t %f \n', 'Weights: ', inh_weights{i}' );
    end
    

%% ======================================================
Final_edges = cell(totalnets,1);
for i = 1 : size(test_edges_act,1)
   activ_edges = test_edges_act{i,1};
   inhibit_edges = test_edges_inh{i,1};
   
   activ_scores = test_scores_act{i,1};
   inhibit_scores = test_scores_inh{i,1};
   
   activ_edges = [activ_edges , activ_scores];
   inhibit_edges = [inhibit_edges , inhibit_scores];
   
   edge = [activ_edges ; inhibit_edges];
   aa = edge(:,[1,2,4]);

   [Ua,~,ix] = unique(aa(:,1:2),'rows');                        % Unique Values & Indices
   L = accumarray(ix, aa(:,3), [], @(x)max(sum(x)~=2*x));       % Eliminate Duplicates (Returns Logical Vector)
   R = accumarray(ix, aa(:,3), [], @max);                       % Maximum Probabilities
   C = [Ua(L,:), R(L)] ;
   [~,iddx] = intersect(aa,C,'rows','stable');
   final_edges = edge(iddx,:);
   
   ts_labels = [ test_labels_act{i,1} ; test_labels_inh{i,1} ] ;
   labels = ts_labels(iddx,:);
   
   Final_edges{i} = [final_edges , labels];                 % Final_edges = [Reg, tar, sign, prob, true_label]
end

load (strcat('genenames_',num2str(num_of_genes),'gene.mat'))
%% ====================== Edge file ==============================
for i = 1 : size(Final_edges,1)
    final_edges = Final_edges{i,1};
    score = final_edges(:,[1,2,4]);
    [~,index] = sort(score(:,3),'descend');

    filename = strcat('all_edges_',num2str(num_of_genes),'_',num2str(i),'.txt');
    fileID = fopen(filename,'w');
    for j = 1:size(index,1)
       idx = index(j,1);
       reg = score(idx,1);
       reg_name = genenames{reg,1};
       target = score(idx,2);
       target_name = genenames{target,1};
       sc = score(idx,3);
       fprintf(fileID,'%s\t%s\t%s\n',reg_name,target_name,num2str(sc));
    end

    fclose(fileID);
end

%% ====================== Probability Matrix =====================
num_of_genes = num_of_genes;
num = num_of_genes*num_of_genes;

for j = 1 : size(Final_edges,1);
    final_edges = Final_edges{j,1};
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

    pr = strcat('prob_', num2str(num_of_genes), '_' , num2str(j),'.mat');
    pr_sign = strcat('prob_sign_', num2str(num_of_genes), '_' , num2str(j),'.mat');
    save(pr, 'probability');
    save(pr_sign, 'prob_sign');
    
end

end
