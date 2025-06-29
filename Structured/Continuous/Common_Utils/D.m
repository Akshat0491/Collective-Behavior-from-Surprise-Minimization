function ret=Dt(A) %derivatibe dot A
%made keeping in mind that A would be a vector of ox1 elements. where o is the order to which generalised coordinates are taken
%returnd D dot A, = D'*A, where D is the derivative matrix
    if nargin<1
        ret=sprintf("D(A) requires a matrix as argument."); %default case, identity matrix
        which Dt
        return; %exit the function

        
    end

    if ~iscell(A)
        temp=zeros(size(A)); 
        for i=1:size(A,1)-1
            temp(i,i+1)=1;
        end
    ret=temp * A;
    else
        B=A;
        for i=1:length(A)-1
            B{i}=B{i+1};
        end
        B{length(B)}={zeros(size(B{length(B)}))};
        ret=B;

    end