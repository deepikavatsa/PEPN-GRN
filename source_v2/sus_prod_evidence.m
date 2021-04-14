function [control,control_sign] = sus_prod_evidence( PreState, op_idx, regulators )

if isempty(regulators)
    total_genes_index = 1:size(PreState,1);
    control_index = total_genes_index';
    control_index(op_idx) = [];
else
    control_index = regulators;
    [~,idx] = ismember(op_idx,regulators);
    if idx ~= 0
        control_index(idx) = [];
    end
end

% Pre allocation
control = zeros(length(control_index),1);
control_sign = zeros(length(control_index),1);

for j = 1 : size(control_index,1)
    con_idx = control_index(j,1);
    con_expr = PreState(con_idx,1);
    if con_expr >= 1
        control(j,1) = con_idx; 
        control_sign(j,1) = 1;
    elseif con_expr == 0
        control(j,1) = con_idx; 
        control_sign(j,1) = 0;
    end
end
    

idx = find(control == 0);
control(idx,:) = [];
control_sign(idx,:) = []; 



end

