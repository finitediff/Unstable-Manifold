clear all; curr_dir = cd; close all; clc;



batch_num = 3;

gamma_vals = linspace(0,2/9,30); 
alpha_vals =  1:1:30;

hold on;
h = xlabel('\gamma');
set(h,'FontSize',22);
h = ylabel('\alpha');
set(h,'FontSize',22);
for gamma = gamma_vals
    p.gamma = gamma;
    for alpha = alpha_vals
        p.alpha = alpha;
        
        file_name = ['batch',num2str(batch_num),'_gray_scott_evans_gamma_', ...
                    num2str(round(p.gamma*1000000)), ...
                    '_alpha_',num2str(round(p.alpha*1000000))];

        try
            ld = retrieve_it(curr_dir,'gray_scott',file_name,'data');
            d = ld.var;
            
            return
            
            if d.wnd1 == 0
                plot(p.gamma,p.alpha,'.r','MarkerSize',18);
            elseif d.wnd1 == 1
                plot(p.gamma,p.alpha,'.k','MarkerSize',18);
            elseif d.wnd2 > 1
                plot(p.gamma,p.alpha,'.g','MarkerSize',18);
                
                cd('..');
                cd('matlab');
                d = evans_batch(d.p,d.s);
            end
            drawnow;
        catch me
            cd(curr_dir);
        end
  
        
    end
end





