function ret=grad(scalar_function,vector) 
    %vector  : array of strings, each been exactly same as the variable it represents, or the direction
    f=sym(scalar_function);
    x=sym(vector);
    gf=gradient(f,x);
    ret=matlabFunction(gf);
end


%example usage


% f=@(x,y) x^2+y^2

% f =

%   function_handle with value:

%     @(x,y)x^2+y^2

% >> ret=grad(f,["x","y"])

% ret =

%   function_handle with value:

%     @(x,y)[x.*2.0;y.*2.0]


% >> ret(1,1)

% ans =

%      2
%      2

% >> ret([1,1])