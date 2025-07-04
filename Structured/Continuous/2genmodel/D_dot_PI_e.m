function ret=D_dot_PI_e(PI,e_mat,o)
    L=size(PI,1)/o;
    b_mat=cell2mat(e_mat);
    Pib=PI*b_mat;
    ret=D(mat2cell(Pib,L*ones(1,o),[1]));
end   