sjcj \frac{a}{b}


Variables:
    o=order to which gnn are truncated, like o=3 if till acceleration
    X= environmental states
    F= how environmental states are evolving
    H= any new environmental aspect that might want to add
    a= action, last entry of the gnn. ex x_tilde={x,xdot,xdotdot,xdotdotdot}, then action would be xdotdotdot
    X_dot = F(X,H,a) + Wx
    VFE:
        order always y,x,v
        y_tilde,mue_tilde_x,mue_tilde_v,PI_tilde_y,PI_tilde_x,PI_tilde_v
    
    stick to updating f(t+1)=f(t)+dt*f'(t)


Y_tilde is a cell containting o cells, each containing n cels, each containing LxT. Suppose o=3, L=5, T=3, N=6






           ___________ __________________________
          |           |                          |
Y_tilde={{};        {{}:            (x1(t=1)  x1(t=2)  x1(t=3))
                                    (x2(t=1)  x2(t=2)  x2(t=3))
                                (   (x3(t=1)  x3(t=2)  x3(t=3))  ) ::(LxT)
                                    (x4(t=1)  x4(t=2)  x4(t=3))
                                    (x5(t=1)  x5(t=2)  x5(t=3))
                     {};
                     {};
                     {};
                     {}::(Nx1)}
         {};   
         {}::(ox1)}







use which(function_name) to fine the location of the function definition


currently implemented where genprocess updates a 2 Dimensional R vector, usually a position vector, and its derivatives. 

