%%
%% state.pp_CH(option)
%%
function pp_CH(obj,option)
    F = dec2bin(obj.F);
    G = dec2bin(obj.G);
    M = dec2bin(obj.M);
    g = dec2bin(obj.g);
    FT = dec2bin(obj.FT);
    GT = dec2bin(obj.GT);
    MT = dec2bin(obj.MT);
    gT = dec2bin(obj.gT);
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
    case 'conj'
        fprintf('// FT (Xp--X)\n');
        disp(FT);
        fprintf('// GT (Zp)\n');
        disp(GT);
        fprintf('// MT (Xp--Z)\n');
        disp(MT);
        fprintf('// gT (Xp--phase)\n');
        disp(gT);
        fprintf('// v (H)\n');
        disp(v);
        fprintf('// |s>\n');
        disp(s);
        fprintf('// global phase\n');
        disp(conj(obj.w)); 
    case 'basis'
        disp(CH2basis(obj));
    end
end
