function ret=update_Y_ext_tilde(Y_tilde,G,R_tilde,noise_params,t)
    %Arguments
    %Y_tilde: Nx1 cell, each cell having ox1 cell, each of which is having LxdlxT
    N=length(Y_tilde);
    o=length(Y_tilde{1});
    L=size(Y_tilde{1}{1},1);
    dL=size(Y_tilde{1}{1},2);
    % parfor i=1:N
    %     for k=1:o
    %     Y_tilde{i}{k}(:,:,t)=G(R_tilde{k}(:,:,t),i,L,dl);
    %     end
    % end

    for k=1:o
        temp=R_tilde{k}(:,:,t);
        parfor i=1:N
            Y_tilde{i}{k}(:,:,t)=G(temp,i,L,dL)
        end
    end

    ret=Y_tilde;
end

% G=@(R_tilde,i,L) sqrt(sum((A(:,1:L)-A(i,1:L)*0).^2, 2));  % returns a column vector of size Lx1, 