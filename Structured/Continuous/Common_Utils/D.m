function ret=D(A) %derivatibe of A/shiftup!
    if isvector(A)
        temp=A;
        o=length(A);
        for i=1:o
            if i==o
                temp(o)=zeros(size(temp(o)));
            else
                temp(i)=temp(i+1);
            end
        end
        ret=temp;
    elseif iscell(A)
        temp=A;
        o=length(A);
        for i=1:o
            if i==o
                temp{o}={zeros(size(temp{o}))};
            else
                temp{i}=temp{i+1};
            end
        end
        ret=temp;
    else
        disp("Unsupported or No input")
    end