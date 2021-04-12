%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'all_pos_edge.mat' and 'all_neg_edge.mat' are the ground truth true positive and 
%% true negative edge matrices. 'all_neg_edge.mat' does not include self edges.
%% 'probability.mat' is the 100x1 column matrix containing probability values of predicted edges 
%% along their respective index. 
%% This script computes the TPR anf FPR in a rankwise manner(i.e., adding each predicted edge at a time)
%% and finally plots ROC and PR curves.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = rankwise_roc_pr_plot(num_of_genes,netnum)


load (strcat('all_pos_edge_',num2str(num_of_genes),'_',num2str(netnum),'.mat'))
load (strcat('all_neg_edge_',num2str(num_of_genes),'_',num2str(netnum),'.mat'))

% load prob matrix
 load (strcat('prob_',num2str(num_of_genes),'_',num2str(netnum),'.mat'))
 load (strcat('prob_sign_',num2str(num_of_genes),'_',num2str(netnum),'.mat'))

% find diagonal indexes
mat = zeros(num_of_genes,num_of_genes);
c = size(mat,1);
diag_idx = 1:c+1:numel(mat);
diag_idx = diag_idx';
clear mat

% sort 'probability' in descending order and delete probability values 0.

[prob,index] = sort(probability,'descend');
prob = prob(~isnan(prob)); % remove NaN entries from prob (which appear at the top of sorted prob matrix)

% remove indexes of self loops from index matrix
for i= 1:size(diag_idx,1)	
	d=diag_idx(i,1);
	index(index==d)=[];
end

prob = prob(prob>0);
length = size(prob,1);
index = index(1:length,1);
sorted_prob_sign = prob_sign(index,:);  % sorted prob_sign according to index vector

% compute TP, FP, TPR, FPR
tp = [0]; fp = [0];
for i = 1:size(prob,1)
   idx = index(i,1) ;
   flag = ismember(idx,all_pos_edge);
   if flag == 1
       tp(1,i+1) = tp(1,i)+1;
       fp(1,i+1) = fp(1,i);
   elseif flag == 0
       tp(1,i+1) = tp(1,i);
       fp(1,i+1) = fp(1,i)+1;
   end
end
[m_pos n_pos]=size(all_pos_edge);
tpr = tp/n_pos;
[m_neg n_neg]=size(all_neg_edge);
fpr = fp/n_neg;


%%=============== GENERATE ROC CURVE ==============================
% appending with max TP and FP values to complete the roc curve at extreme points

TP = [tp,n_pos];
FP = [fp,n_neg];
TPR = [tpr,1];
FPR = [fpr,1];

AUROC=abs(trapz(FPR,TPR)); 	% Here, first argument should be the 'x' variable of curve i.e., FPR
AUROC = (round(AUROC,2));

h2=figure('Position', [300 300 390 200]);
if all(TPR(:)==1) && all(FPR(:)==0)
    plot(FPR,TPR,'*')
else
    plot(FPR,TPR,'o-','Linewidth',2,'markersize',5,'markerfacecolor','g')
end
set(gca,'FontSize',13);
xlabel('FPR','FontSize', 20);
ylabel('TPR','FontSize', 20);
titlename='PEPN:ROC plot';
title(titlename,'Interpreter','none','FontSize',20);
leg2=legend(sprintf('AUROC %0.2f',AUROC));
set(leg2,'FontSize',15,'Location','southeast');
set(h2, 'PaperPositionMode','auto')
print(strcat('rocplot_',num2str(num_of_genes),'_',num2str(netnum)),'-depsc')


AUROC = num2str(AUROC);
AUROC

%=====================   GENERATE PR PLOT   =================================
TP = tp;
FP = fp;
TPR = tpr;
FPR = fpr;

for i=1:size(TP,2)
   if (TP(i)+FP(i))==0
       Precision(i) = 1;
   else
       Precision(i) = TP(i)/(TP(i)+FP(i));
   end
   Recall(i)=TPR(i);
   Fscore(i)=2*( (Precision(i)*Recall(i)) / (Precision(i)+Recall(i)) );
end

AUPR=abs(trapz(Recall,Precision)); % note that first write 'x' variable of curve i.e., Recall here and then 'y' variable to get correct result
AUPR = (round(AUPR,2));

h1=figure('Position', [300 300 390 200]);
plot(Recall,Precision,'o-','Linewidth',2,'markersize',5,'markerfacecolor','g');
set(gca,'Fontsize',13);
xlabel('Recall','FontSize', 20);
ylabel('Precision','FontSize', 20);
titlename='PEPN:PR plot';
title(titlename,'Interpreter','none','FontSize',20);
leg1=legend(sprintf('AUPR %0.2f',AUPR));
set(leg1,'FontSize',15,'Location','northeast');
set(h1, 'PaperPositionMode','auto')
print(strcat('prplot_',num2str(num_of_genes),'_',num2str(netnum)),'-depsc')

AUPR = num2str(AUPR);   % done to remove trailing zeros from the AUPR
AUPR

end
