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
            obj.set_FT(i,i,1);
            obj.set_GT(i,i,1);
        end
        obj.w = 1;
    case 'rand' % generate by applying random gates!
        obj.CH_init('zero');
        obj.s = uint16(randi(obj.len,1,1) - 1);
        for i = 1:30
            obj.CH_gate('rand',[-1,-1]);
        end
    otherwise
        fprintf('error CH_init: invalid option string.\n');
    end
end
