sjcj \frac{a}{b}


Variables:
    X= environmental states
    F= how environmental states are evolving
    H= any new environmental aspect that might want to add
    a= action, last entry of the gnn. ex x_tilde={x,xdot,xdotdot,xdotdotdot}, then action would be xdotdotdot
    X_dot = F(X,H,a) + Wx
    VFE:
        order always y,x,v
        y_tilde,mue_tilde_x,mue_tilde_v,PI_tilde_y,PI_tilde_x,PI_tilde_v
    
    stick to updating f(t+1)=f(t)+dt*f'(t)



use which(function_name) to fine the location of the function definition
