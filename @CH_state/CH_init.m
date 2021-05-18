%%
%%  state.CH_init(option)
%%   
function CH_init(obj,option)
    % todo: add rand option
    % todo: combine init and reset? maybe inefficient and redundant
    switch option
    case 'zero' % |0^n> 
        for i = 1:obj.len
            obj.set_F(i,i,1);
            obj.set_G(i,i,1);
        end
        obj.w = 1;
    otherwise
        fprintf('error CH_init: invalid option string.\n');
    end
end


