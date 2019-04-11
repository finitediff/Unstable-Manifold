function U_n = deval_loaded(dom, s,p)
    if (isnumeric(dom) == false)
        temp = dom;
        dom = zeros(1,1);
        dom(1) = temp;
    end
        
    U_n = zeros(4,length(dom));
    for j = 1:length(dom)
        if dom(j) >= s.L
            if dom(j) <= s.R
                temp = soln(dom(j),s);
            else
                temp = soln(s.R,s);
            end
        else
            temp = soln(s.L,s);
        end
       U_n(:,j) = [1-temp(1);temp(1);temp(2);temp(3)];
    end

end

