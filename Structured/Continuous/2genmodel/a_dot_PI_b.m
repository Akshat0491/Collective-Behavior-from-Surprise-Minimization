function ret=a_dot_PI_b(a,PI,b)
    % a       : cell(o,1), each cell having Lx1 vector
    % PI      : cell(1,dl), each cell having oLxoL matrix
    % b       : cell(o,1), each cell having Lxdl matrix
    % dL=length(PI);  %will think of this later, it is set to 1 for now...
    
    o=size(a,1);
    L=size(a{1},1);
    % ret=cell(o,1);
    b_mat=cell2mat(b);
    Pib=b_mat'*PI*b_mat;
    ret=Pib;
    % b_cell=mat2cell(Pib,L*ones(1,o),[1]);
    % for i=1:o
    %     ret{i}=a{i} .* b_cell{i};
    % end

end   