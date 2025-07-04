function ret=D(A) %derivatibe of A/shiftup!
    if ~iscell(A)
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
    elseif iscell(A) %made mainly for mu_tilde
        temp=A;
        o=length(A);
        for i=1:o
            if i==o
                tempo=zeros(size(A{1}));

                temp{o}=tempo;
            else
                temp{i}=temp{i+1};
            end
        end
        ret=temp;
    else
        disp("Unsupported or No input")
    end