%%
%%  stab_new = stab.CH_gate(g,q)
%%                                             
function stab_new = CH_gate(obj,g,param)
    % todo: makes shallow copy for now and mutates original state; 
    % todo: make class methods void functions: mutation->void ; no mutation->return new resulting state
    % may want to deepcopy depending on the search algorithm (shallowcopy sufficient for SA)
    switch g
    case 'SL'
        stab_new = CH_SL(param,obj);
    case 'SR'
        stab_new = CH_SR(param,obj);
    case 'CZL'
        stab_new = CH_CZL(param,obj);
    case 'CZR'
        stab_new = CH_CZR(param,obj);
    case 'CXL'
        stab_new = CH_CXL(param,obj);
    case 'CXR'
        stab_new = CH_CXR(param,obj);
    case 'HL'
        stab_new = CH_HL(param,obj);
    case 'STL' % S transpose used for updating (Uc)^T
        stab_new = CH_HL(param,obj);
    case 'STR'
        stab_new = CH_HL(param,obj);
    otherwise
        fprintf('error CH_gate: invalid gate name.\n');
    end
end

% apply C-type gates from left
function stab_new = CH_SL(param,stab)
    q = param(1);
    stab.M(q,1) = bitxor(stab.M(q,1),stab.G(q,1));
    %for p = 1:stab.len
    %    stab.set_M(p,q,bitxor(stab.get_M(p,q),stab.get_G(p,q)));
    %end
    new_gq = stab.get_g(q);
    new_gq = bitset(new_gq,2,~bitxor(bitget(new_gq,1),bitget(new_gq,2)));
    new_gq = bitset(new_gq,1,~bitget(new_gq,1));
    stab.set_g(q,new_gq); 
    stab_new = stab;
end

function stab_new = CH_CZL(param,stab)
    q = param(1);
    r = param(2);
    %for p = 1:stab.len
    %    stab.set_M(q,p,bitxor(stab.get_M(q,p),stab.get_G(r,p)));
    %    stab.set_M(r,p,bitxor(stab.get_M(r,p),stab.get_G(q,p)));
    %end
    stab.M(q,1) = bitxor(stab.M(q,1),stab.G(r,1));
    stab.M(r,1) = bitxor(stab.M(r,1),stab.G(q,1));
    stab_new = stab;
end

function stab_new = CH_CXL(param,stab)
    q = param(1);
    r = param(2);
    %fprintf('g_q:%d, g_r:%d\n, g_q+g_r:%d\n',stab.get_g(q),stab.get_g(r),mod(stab.get_g(q)+stab.get_g(r),4));
    stab.set_g(q,mod(stab.get_g(q)+stab.get_g(r),4));
    %fprintf('new g_q: %d\n',stab.get_g(q));
    % todo: optimize bit ops
    if parity(bitand(stab.M(q),stab.F(r)))
        new_gq = stab.get_g(q);
        new_gq = bitset(new_gq,2,~bitget(new_gq,2));
        stab.set_g(q,new_gq); 
    end
    %fprintf('M_q: %s,F_r: %s\n',dec2bin(stab.M(q)),dec2bin(stab.F(r)));
    %fprintf('par: %d\n',parity(bitand(stab.M(q),stab.F(r))));
    %fprintf('new new g_q: %d\n',stab.get_g(q));
    %for p = 1:stab.len
    %    stab.set_G(r,p,bitxor(stab.get_G(r,p),stab.get_G(q,p)));
    %    stab.set_F(q,p,bitxor(stab.get_F(q,p),stab.get_F(r,p)));
    %    stab.set_M(q,p,bitxor(stab.get_M(q,p),stab.get_M(r,p)));
    %end
    stab.G(r,1) = bitxor(stab.G(r,1),stab.G(q,1));
    stab.F(q,1) = bitxor(stab.F(q,1),stab.F(r,1));
    stab.M(q,1) = bitxor(stab.M(q,1),stab.M(r,1));
    stab_new = stab;
end

% apply C-type gates from right 
function stab_new = CH_SR(param,stab)
    q = param(1);
    for p = 1:stab.len
        stab.set_M(p,q,bitxor(stab.get_M(p,q),stab.get_F(p,q)));
        % todo: optimize bit ops
        if stab.get_F(p,q)
            new_gp = stab.get_g(p);
            new_gp = bitset(stab.get_g(p),2,~bitxor(bitget(new_gp,1),bitget(new_gp,2)));
            new_gp = bitset(new_gp,1,~bitget(new_gp,1));
            stab.set_g(p,new_gp); 
        end
    end
    stab_new = stab;
end

function stab_new = CH_CZR(param,stab)
    q = param(1);
    r = param(2);
    for p = 1:stab.len
        stab.set_M(p,q,bitxor(stab.get_M(p,q),stab.get_F(p,r)));
        stab.set_M(p,r,bitxor(stab.get_M(p,r),stab.get_F(p,q)));
        % todo: optimize bit ops
        if stab.get_F(p,q) && stab.get_F(p,r)
            new_gp = stab.get_g(p);
            new_gp = bitset(new_gp,2,~bitget(new_gp,2));
            stab.set_g(p,new_gp); 
        end
    end
    stab_new = stab;
end

function stab_new = CH_CXR(param,stab)
    q = param(1);
    r = param(2);
    for p = 1:stab.len
        stab.set_G(p,q,bitxor(stab.get_G(p,q),stab.get_G(p,r)));
        stab.set_F(p,r,bitxor(stab.get_F(p,r),stab.get_F(p,q)));
        stab.set_M(p,q,bitxor(stab.get_M(p,q),stab.get_M(p,r)));
    end
    stab_new = stab;
end

% apply H gate on the left (see prop'n 4 from "Simulation of qc by low-rank stab decomp" paper)
function stab_new = CH_HL(param,stab)
    p = param(1);
    % step 1: move HLq inside UcUh layer; calculate  constants
    t = bitxor(stab.s,bitand(stab.G(p,1),stab.v));
    %fprintf('t:');disp(dec2bin(t));
    u = bitxor(bitxor(stab.s,bitand(stab.M(p,1),stab.v)),bitand(stab.F(p,1),bitcmp(stab.v)));
    %fprintf('u:');disp(dec2bin(u));
    a = parity(bitand(bitand(stab.G(p,1),bitcmp(stab.v)),stab.s));
    %fprintf('a:');disp(a);
    b = parity( ...
            bitxor( ...
                bitand(bitand(stab.M(p,1),bitcmp(stab.v)),stab.s), ...
                bitand(bitand(stab.F(p,1),stab.v),bitxor(stab.M(p,1),stab.s)) ...
            ) ...
        );
    %fprintf('b:');disp(b);
    d = mod(stab.get_g(p)+a+a+b+b,4);
    %fprintf('d:');disp(d);

    if t ~= u
        superpos2circuit(t,u,a,d,stab);
    else
        stab.w = stab.w * 2.^(-0.5) * ((-1).^double(a) + ((1i).^double(stab.get_g(p))) * ((-1).^double(b)));
        stab.s = t;
    end
    stab_new = stab;
end







