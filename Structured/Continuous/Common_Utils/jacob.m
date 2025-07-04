function ret=jacob(vector_function,vector) 
    %vector  : array of strings, each been exactly same as the variable it represents, or the direction
    f=sym(vector_function);
    x=sym(vector);
    jacobi=jacobian(f,x);
    sz=size(jacobi);
    ret=cell(sz);
    for i=1:sz(1)
        for j=1:sz(2)
            func=matlabFunction(jacobi(i,j));
            if nargin(func)==0
                func1=@(a) func();
            % elseif nargin(func)==1
            %     func1=@(a,b) func(a);
            else
                func1=func;
            end
            ret{i,j}=func1;

            % ret{i,j}=func;
        end
    end

    % ret=jacob;
end


















% Utility to wrap a function handle with zero arguments to accept two arguments
function wrapped = wrapZeroArgFunc(func)
    wrapped = @(a,b) func();
end




















