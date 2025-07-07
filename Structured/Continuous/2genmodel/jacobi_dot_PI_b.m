function ret=jacobi_dot_PI_b(jacobi,PI,b)
    % jacobi  : cell(o,o), each cell having Lx1 vector
    % PI      : cell(1,dl), each cell having oLxoL matrix
    % b       : cell(o,1), each cell having Lxdl matrix
    % dL=length(PI);  %will think of this later, it is set to 1 for now...
    
    sz=size(jacobi);
    L=size(b{1},1);
    o=sz(1);
    d=sz(2);
    ret=cell(d,1);
    b_mat=cell2mat(b);
    Pib=PI*b_mat;
    b_cell=mat2cell(Pib,L*ones(1,o),[1]);
    for i=1:d
        ret{i}=zeros(L,1);
        
        for j=1:o
            ret{i}=ret{i} + jacobi{j,i} .* b_cell{j};
        end
    end

end   