a=sym("a","real");
mx=sym("mx","real");
mv=sym("mv","real");
ma=sym("ma","real");
yx=sym("yx","real");
yv=sym("yv","real");
eta=sym("eta","real");
x2=sym("x2","real");
x1=sym("x1","real");
v2=sym("v2","real");
v1=sym("v1","real");

ey=[x2-x1-mx;v2-v1-mv]
ex=[mv-(-a*(mx-eta));ma+(a*mv); a*ma]
PI_x=init_PI_tilde(3,1,[1],[1]);
PI_y=init_PI_tilde(2,1,[1],[1]);

F=dot(ey,PI_y*ey) + dot(ex,PI_x*ex)
gF_v=gradient(F,[v1]);
gradfunc=matlabFunction(gF_v)
% gF_act=