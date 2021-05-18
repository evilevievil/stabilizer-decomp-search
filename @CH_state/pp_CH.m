%%
%% state.pp_CH(option)
%%
function pp_CH(obj,option)
    % todo: helper to remove duplicate code for F,G,M
    F = dec2bin(obj.F);
    G = dec2bin(obj.G);
    M = dec2bin(obj.M);
    g = dec2bin(obj.g);
    v = dec2bin(obj.v);
    s = dec2bin(obj.s);
    switch option
    case 'ch'
        fprintf('// F (Xp--X)\n');
        disp(F);
        fprintf('// G (Zp)\n');
        disp(G);
        fprintf('// M (Xp--Z)\n');
        disp(M);
        fprintf('// g (Xp--phase)\n');
        disp(g);
        fprintf('// v (H)\n');
        disp(v);
        fprintf('// |s>\n');
        disp(s);
        fprintf('// global phase\n');
        disp(obj.w); 
    end
end
