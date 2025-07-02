function ret=update_sense(Y_ext_tilde,G,R_tilde,noise_params,PI_y,t)
    Y_ext_tilde=update_Y_ext_tilde(Y_ext_tilde,G,R_tilde,noise_params,t);
    PI_y       =update_PI_tilde(PI_y);



    vars = ["Y_ext_tilde", "PI_y"];
    vals = [{Y_ext_tilde}, {PI_y}];
    ret  = dictionary(vars,vals);
end