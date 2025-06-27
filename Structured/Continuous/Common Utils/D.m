function ret=Dt(A) %derivatibe dot A
%made keeping in mind that A would be a vector of ox1 elements. where o is the order to which generalised coordinates are taken
%returnd D dot A, = D'*A, where D is the derivative matrix
    if nargin<1
        ret=sprintf("D(A) requires a matrix as argument."); %default case, identity matrix
        which Dt
        return; %exit the function

        
    end
    temp=zeros(size(A)); 
    for i=1:size(A,1)-1
        temp(i,i+1)=1;
    end
ret=temp * A;
end