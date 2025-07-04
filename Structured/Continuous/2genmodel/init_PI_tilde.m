function ret=init_PI_tilde(o,L,gamma_vecdL,lambda_vecdL)
    %return a cell(1,dl), each cell containing olxol matrix. 

    %gamma and lambda _vecdl are 1xdl vectors
    dL=length(gamma_vecdL); %by default dL=1
    PI=cell(1,dL);

    parfor j=1:dL
        Spatial_term_j=inv(gamma_vecdL(j)*eye(L));
        Temporal_term_j=spm_DEM_R(o,lambda_vecdL(j));
        % PI{j}=kron(Spatial_term_j,Temporal_term_j);
        PI{j}=kron(Temporal_term_j,Spatial_term_j); %kept like this intentionally, w.r.t how I just flatten out mu_tilde later
    end

    ret=PI{1};
end




%" WeparameterizetheL×LspatialcovarianceΣωthroughitsprecisionmatrixΠω,asadiagonalmatrixwhoseentriesare 284 givenbyasingleprecision(inversevariance)Γω similar to heins et al2024" 