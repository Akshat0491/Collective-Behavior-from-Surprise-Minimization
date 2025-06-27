function ret=vfe(y_tilde,mue_tilde_x,mue_tilde_v,PI_tilde_y,PI_tilde_x,PI_tilde_v) %add methods for gradients w.r.t diff parameters
% VFE - Variational Free Energy

epsilon_tilde_y = y_tilde - g_tilde(mue_tilde_x,mue_tilde_v);
epsilon_tilde_x = D(mue_tilde_x)-f_tilde(mue_tilde_x,mue_tilde_y);
epsilon_tilde_v =mue_tilde_v-eta_tilde;

ret=0.5 * dot(epsilon_tilde_y,PI_tilde_y*epsilon_tilde_y) + ...
    0.5 * dot(epsilon_tilde_x,PI_tilde_x*epsilon_tilde_x) + ...
    0.5 * dot(epsilon_tilde_v,PI_tilde_v*epsilon_tilde_v);

end