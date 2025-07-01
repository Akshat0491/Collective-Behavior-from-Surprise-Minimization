function ret=a_dot_PI_b(a,PI,b,t)
    % a  : cell(o,1), each cell having LxdlxT matrix
    % PI : cell(1,dl), each cell having oLxoL matrix
    % b  : cell(o,1), each cell having LxdlxT matrix
    dL=length(PI);
    a_mat=a;
    if iscell(a)
    a_mat=cell2mat(a);
    else
        disp("a is not a cell")
    end

    b_mat=b;
    if iscell(b)
    b_mat=cell2mat(b);
    else
        disp("b is not a cell")
    end
    prod=cell(1,dL);
    parfor j=1:dL 
        prod{j}=a_mat(:,j,t)'*PI{j}*b_mat(:,j,t)
    end
    ret=cell2mat(prod);
end