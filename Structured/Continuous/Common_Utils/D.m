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
                tempo=cell(length(A{1}),1);
                temptemp=temp{1};
                parfor l=1:length(A{1})
                    tempo{l}=zeros(size(temptemp));
                end
                temp{o}=tempo;
            else
                temp{i}=temp{i+1};
            end
        end
        ret=temp;
    else
        disp("Unsupported or No input")
    end