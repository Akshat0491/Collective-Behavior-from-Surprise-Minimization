function ret=shift_down_PI_e(PI,e)
    function ret=jacobi_dot_PI_b(jacobi,PI,b_mat,t)
    % a  : cell(o,1), each cell having LxdlxT matrix
    % PI : cell(1,dl), each cell having oLxoL matrix
    % b  : cell(o,1), each cell having LxdlxT matrix
    dL=length(PI);
    prod=cell(1,dL);
    parfor j=1:dL 
        prod{j}=PI{j}*b_mat(:,j,t)
    end
    ret=D(prod);
end   