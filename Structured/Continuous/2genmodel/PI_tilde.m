%write a function that inputs: gamma,lambda,L
function ret=PI_tilde(L,gamma,lambda)
    Spatial_term=inv(gamma*eye(L));
    Temporal_term=spm_DEM_R(L,lambda);
    ret=kronecker(Spatial_term,Temporal_term);
end




%" WeparameterizetheL×LspatialcovarianceΣωthroughitsprecisionmatrixΠω,asadiagonalmatrixwhoseentriesare 284 givenbyasingleprecision(inversevariance)Γω similar to heins et al2024" 