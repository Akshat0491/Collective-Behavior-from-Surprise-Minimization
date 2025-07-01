function ret=init_PI_tilde(o,L,gamma_vecdL,lambda_vecdL)
    %return a cell(1,dl), each cell containing olxol matrix. 

    %gamma and lambda _vecdl are 1xdl vectors
    dL=length(gamma_vecdL);
    PI=cell(1,dL);

    parfor j=1:dL
        Spatial_term_j=inv(gamma_vecdL(j)*eye(L));
        Temporal_term_j=spm_DEM_R(o,lambda_vecdL(j));
        PI{j}=kron(Spatial_term_j,Temporal_term_j);
    end

    ret=PI;
end




%" WeparameterizetheL×LspatialcovarianceΣωthroughitsprecisionmatrixΠω,asadiagonalmatrixwhoseentriesare 284 givenbyasingleprecision(inversevariance)Γω similar to heins et al2024" 