%%
%%  state.CH_init(option)
%%   
function CH_init(obj,option)
   % todo: add rand option
   % all zero state
   row = uint8(1);
   for i = 1:8
      set(obj,'F',i,i,1);
      set(obj,'G',i,i,1);
   end
   set(obj,'w',1);
end


