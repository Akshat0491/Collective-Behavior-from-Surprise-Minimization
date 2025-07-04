function ret=get_curr_sense(Y_ext_tilde,t)
    N=length(Y_ext_tilde);
    o=length(Y_ext_tilde{1});
    ret=cell(N,1);
    parfor i=1:N
        ret{i}=cell(o,1);
        for k=1:o
            ret{i}{k}=Y_ext_tilde{i}{k}(:,:,t);
        end
    end
end
