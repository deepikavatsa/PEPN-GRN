function [SMOTE_samples, SMOTE_labels] = smote(samples, labels)

idx = labels == 1;                      % Logical index
min_samples = samples(idx,:);           % minority class samples

idx = labels == 0;
maj_samples = samples(idx,:);           % majority class samples

maj_num = size(maj_samples, 1);         % num of majority samples
min_num = size(min_samples, 1);         % num of minority samples
num_samples_to_gen = maj_num - min_num; % total num of synthetic sample to generate

%---------------------------------------------------------
% SMOTE Algorithm
% num_samples_to_gen / min_num = this value will be in fraction
% N = floor(num_samples_to_gen / min_num) --> N will be the floor value of
% that fraction. Therefore, N * min_num will produce less samples than
% needed. rest_samples denote num of more samples to produce

N = floor(num_samples_to_gen / min_num);        % N = num of synthetic samples to gen from each minority sample
num = N * min_num;                              % how many will be generated?
rest_samples = num_samples_to_gen - num;        % rest samples to be generated such that total num of syn samples are generated

k = 5;                                          % number of nearest neighbours

numfeatures = size(min_samples,2);              % number of features
syn_samples = zeros(num_samples_to_gen , numfeatures);

% I = k nearest neighbours indexes for each minority class sample
I = nearestneighbour(min_samples', min_samples', 'NumberOfNeighbours', k+1);
I = I';

j = 1;
for i = 1 : min_num
    nnarray = I(i,2:end);
    N1 = N;
    while N1 ~= 0
    nnidx = randi([1 k] , 1);
    nn = nnarray(1 , nnidx);                    % selecting one nearest neighbour at random
    
        for f = 1 : numfeatures                 % generate a synthetic sample
            diff = min_samples(nn , f) - min_samples(i , f);
            gap = rand;
            syn_samples(j,f) = min_samples(i , f) + (gap * diff);
        end
        j = j+1;
        N1 = N1 - 1;
    end
    
end

for i = 1 : rest_samples
    nnarray = I(i,2:end);
    nnidx = randi([1 k] , 1);
    nn = nnarray(1 , nnidx);
    
        for f = 1 : numfeatures
            diff = min_samples(nn , f) - min_samples(i , f);
            gap = rand;
            syn_samples(j,f) = min_samples(i , f) + (gap * diff);
        end
        j = j+1;
end


nrow = size(syn_samples,1);
syn_label = ones(nrow, 1);
    
SMOTE_samples = [samples ; syn_samples];
SMOTE_labels = [labels ; syn_label];

end
    
