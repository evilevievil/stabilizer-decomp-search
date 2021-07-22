%%
%%  state.CH_init(option)
%%   
function CH_init(obj,option)
    switch option
    case 'zero' % |0^n> 
        for i = 1:obj.len
            obj.set_F(i,i,1);
            obj.set_G(i,i,1);
            obj.set_FT(i,i,1);
            obj.set_GT(i,i,1);
        end
        obj.w = 1;
    case 'rand' % generate random state using random circuit
        obj.CH_init('zero');
        bit_max = 2.^(obj.len);
        obj.s = cast(randi(bit_max,1,1) - 1,'like',const.init_uint);
        for i = 1:100
            obj.CH_gate('rand',[-1,-1]);
        end
    otherwise
        fprintf('error CH_init: invalid option string.\n');
    end
end
