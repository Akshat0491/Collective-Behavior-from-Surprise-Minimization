function ret=make_jacobi_matrix_action(i,t,o,L,jacobian_func_int_tilde,vals,RR_tilde)
    %vals : o*L x1 array. L;L;L...
            vals=reshape(vals,[L,o])';   %makes in oxL
            sz=size(jacobian_func_int_tilde);
        jacobi_mat=cell(sz);    %each cell shal get Lx1 thing
        for m=1:sz(1)
            for n=1:sz(2)
                m,n
                tempp=jacobian_func_int_tilde{m,n}(RR_tilde{m}(:,1,t-1)+eps,RR_tilde{m}(:,2,t-1))
                if isscalar(tempp)
                    jacobi_mat{m,n}=tempp*ones(1,L)';
                else
        jacobi_mat{m,n}=tempp';
                end
            end
        end
        ret=jacobi_mat;
    end
    














