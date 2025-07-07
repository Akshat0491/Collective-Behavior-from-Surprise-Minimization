function ret=make_jacobi_matrix(o,L,jacobian_func_int_tilde,vals)
    %vals : o*L x1 array. L;L;L...
            vals=reshape(vals,[L,o])';   %makes in oxL
            sz=size(jacobian_func_int_tilde);
        jacobi_mat=cell(sz);    %each cell shal get Lx1 thing
        for m=1:sz(1)
            for n=1:sz(2)
                if isscalar(jacobian_func_int_tilde{m,n}(vals(m,:)))
                    jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(vals(m,:))*ones(1,L)';
                else
        jacobi_mat{m,n}=jacobian_func_int_tilde{m,n}(vals(m,:))';
                end
            end
        end
        ret=jacobi_mat;
end
