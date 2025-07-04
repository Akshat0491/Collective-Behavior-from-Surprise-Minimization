function ret=make_jacobi_matrix(o,L,jacobian_func_int_tilde,vals)
    %vals : o*L x1 array. L;L;L...
            vals=reshape(vals,[L,o])';   %makes in oxL
        jacobi_mat=cell(o,o);    %each cell shal get Lx1 thing
        for m=1:o
            for n=1:o
                if isscalar(jacobian_func_int_tilde{m,n}(vals(m,:)))
                    jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(vals(m,:))*ones(1,L)';
                else
        jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(vals(m,:))';
                end
            end
        end
        ret=jacobi_mat;
    end























    