% a=sym("a","real");
% mx=sym("mx","real");
% mv=sym("mv","real");
% ma=sym("ma","real");
% yx=sym("yx","real");
% yv=sym("yv","real");
% eta=sym("eta","real");
% x2=sym("x2","real");
% x1=sym("x1","real");
% v2=sym("v2","real");
% v1=sym("v1","real");
% az1=sym("az1","real");
% a1=sym("a","real");
% z1=sym("z1","real");
% z2=sym("z2","real");
% zv1=sym("zv1","real");
% mu_theta=sym("mu_theta","real");




% ey=[atan2d(zv1,v1)-mx; 
%     ((v1*az1 - zv1*a1)/(az1^2 + a1^2))-mv]

% % ey=[yx-mx; yv-mv]

% ex=[mv-(-a*(mx-eta));ma+(a*mv); a*ma]

% PI_x=init_PI_tilde(3,1,[1],[1]);

% PI_y=init_PI_tilde(2,1,[1],[1]);

% F=dot(ey,PI_y*ey) + dot(ex,PI_x*ex)
% gF_v=gradient(F,[mx,mv,ma])
% gradfunc=matlabFunction(gF_v)
% % gF_act=
% % latex(gF_v)

% theta=sym("theta","real");
% mu_theta=sym("mu_theta","real");
% theta_dot=sym("theta_dot","real");
% mu_theta_dot=sym("mu_theta_dot","real");
alph=sym("alph","real");
eta_x=sym("eta_x","real");
v0x=sym("v0x","real");
v0y=sym("v0y","real");
vx=sym("vx","real");
vy=sym("vy","real");
% mu_theta_dot_dot=sym("mu_theta","real");
% % mu_theta=sym("mu_theta","real");

% y_phi=sym("y_phi","real");
mu_phi=sym("mu_phi","real");
mu_phi_dot=sym("mu_phi_dot","real");

y_phi=vx*v0x + vy*v0y;

ey=[y_phi     - mu_phi];

ex=[(mu_phi_dot+(alph*(mu_phi-eta_x)));(-alph*(mu_phi_dot))]

noise_params=[1,1,1];
PI_x=diag(noise_params(1:2));
PI_y=noise_params(3);

f=dot(ey,PI_y*ey) + dot(ex,PI_x*ex)
gFm=gradient(f,[vx,vy])
ret=matlabFunction(gFm)
% latex(gFm)










