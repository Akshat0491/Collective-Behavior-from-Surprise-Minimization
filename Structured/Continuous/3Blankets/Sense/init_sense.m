function ret=init_sense(N,o,L,dL,T,gamma_vecdL_y,lambda_vecdL_y)
    Y_ext_tilde =init_Y_tilde(N,o,L,dL,T);
    PI_y        =init_PI_tilde(o,L,gamma_vecdL_y,lambda_vecdL_y);



    vars = ["Y_ext_tilde", "PI_y"];
    vals = [{Y_ext_tilde}, {PI_y}];
    ret  = dictionary(vars,vals);


end