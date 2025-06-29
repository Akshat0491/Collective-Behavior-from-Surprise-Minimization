function ret=Dt(A) %derivatibe dot A/shiftdown!
    if isvector(A)
        temp=A;
        o=length(A);
        for i=1:o
            if i==1
                temp(i)=zeros(size(A(1)));
            else
                temp(i)=A(i-1);
            end
        end
        ret=temp;
    elseif iscell(A)
        temp=A;
        o=length(A);
        for i=1:o
            if i==1
                temp{i}={zeros(size(A{1}))};
            else
                temp{i}=temp{i+1};
            end
        end
        ret=temp;
    else
        disp("Unsupported or No input")
    end