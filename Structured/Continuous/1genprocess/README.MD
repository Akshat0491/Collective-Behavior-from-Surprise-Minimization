currently very shallow

pseudocode for genprocess_gnn

x_tilde = { [NxT];
            [NxT];
            [NxT];  }
o=length(t) = 3

F_tilde =    {  [zeros(N,T)];
                [zeros(N,T)]
initialised as  [zeros(N,T)];  }

F_tilde{o}=F(x_tilde,action)-->
          { [zeros(N,T)];
            [zeros(N,T)];    
            [f(x,u):Nx1];  }

i=1
x_tilde_temp=x_tilde{o}=[acc:NxT]
x_tilde_temp(:,t+1)=x_tilde_temp(:,t)+ dt * F_tilde(o-i+1);

[Nx(T=t+1)]=[Nx(T=t)] + dt .* F_tilde(o)
...




Future avenue: seperate out f_tilde generation. currentl F_tilde generation is part of the genprocess_gnn function



\frac{a}{b} 




                                Flow of Functions
      

                              
      ___________________          
      | init_genprocess |          
--------------------------------   
|                             |    
| init_R_tilde & init_F_tilde |   
|                             |    
-------------------------------   

              |
              |
              |
              |
              |
              Y
      ___________________          
      | update_genprocess |          
--------------------------------   
|                             |    
| update_F_tilde  (action)    |
| update_R_tilde  (F)         |
|                             |    
-------------------------------   
