%%
%%  stab_new = stab.CH_gate(g,q)
%%                                             
function stab_new = CH_gate(obj,g,param)
    % todo: makes shallow copy for now and mutates original state; 
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
    otherwise
        fprintf('error CH_gate: invalid gate name.\n');
    end
end

% apply C-type gates from left
function stab_new = CH_SL(param,stab)
    q = param(1);
    for p = 1:stab.len
        stab.set_M(p,q,bitxor(stab.get_M(p,q),stab.get_G(p,q)));
    end
    new_gq = stab.get_g(q);
    new_gq = bitset(new_gq,2,~bitxor(bitget(new_gq,1),bitget(new_gq,2)));
    new_gq = bitset(new_gq,1,~bitget(new_gq,1));
    stab.set_g(q,new_gq); 
    stab_new = stab;
end

function stab_new = CH_CZL(param,stab)
    q = param(1);
    r = param(2);
    for p = 1:stab.len
        stab.set_M(q,p,bitxor(stab.get_M(q,p),stab.get_G(r,p)));
        stab.set_M(r,p,bitxor(stab.get_M(r,p),stab.get_G(q,p)));
    end
    stab_new = stab;
end

function stab_new = CH_CXL(param,stab)
    q = param(1);
    r = param(2);
    stab.set_g(q,mod(stab.get_g(q)+stab.get_g(r),4));
    % todo: optimize bit ops
    if parity(bitand(stab.M(q),stab.F(r)))
        new_gq = stab.get_g(q);
        new_gq = bitset(new_gq,2,~bitget(new_gq,2));
        stab.set_g(q,new_gq); 
    end
    for p = 1:stab.len
        stab.set_G(r,p,bitxor(stab.get_G(r,p),stab.get_G(q,p)));
        stab.set_F(q,p,bitxor(stab.get_F(q,p),stab.get_F(r,p)));
        stab.set_M(q,p,bitxor(stab.get_M(q,p),stab.get_M(r,p)));
    end
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
    a = uint8(parity(bitand(bitand(stab.G(p,1),bitcmp(stab.v)),stab.s)));
    %fprintf('a:');disp(a);
    b = uint8(parity( ...
            bitxor( ...
                bitand(bitand(stab.M(p,1),bitcmp(stab.v)),stab.s), ...
                bitand(bitand(stab.F(p,1),stab.v),bitxor(stab.M(p,1),stab.s)) ...
            ) ...
        ));
    %fprintf('b:');disp(b);
    d = mod(stab.get_g(p)+a+a+b+b,4);
    %fprintf('d:');disp(d);

    if t ~= u
    % step 2: compute strings y,z (differ by one bit q), string v0,v1 representing Vc
        t_oxr_u = bitxor(t,u);
        %fprintf('t_oxr_u:');disp(dec2bin(t_oxr_u));
        v0 = bitand(t_oxr_u,bitcmp(stab.v));
        %fprintf('v0:');disp(dec2bin(v0));
        v1 = bitand(t_oxr_u,stab.v);
        %fprintf('v1:');disp(dec2bin(v1));
        q = 1;
        while ~bitget(t_oxr_u,q)
            q = q+1;
        end
        if bitget(t,q)
            y = bitxor(u,bitset(uint8(0),q));
            %fprintf('y:');disp(dec2bin(y));
        else
            y = t;
            %fprintf('y:');disp(dec2bin(y));
        end
        
        % update Uc = UcVc
        if bitsum(t_oxr_u) > 1
            if v0 
                for pos = 1:stab.len
                    if bitget(v0,pos) && (pos ~= q)
                        %fprintf('curr v0 pos:');disp(pos);
                        stab.CH_gate('CXR',[q,pos]);
                    end
                end
                for pos = 1:stab.len
                    if bitget(v1,pos)
                        %fprintf('curr v1 pos:');disp(pos);
                        stab.CH_gate('CZR',[q,pos]);
                    end
                end
            elseif v1
                for pos = 1:stab.len
                    if bitget(v1,pos) && (pos ~= q)
                        stab.CH_gate('CXR',[pos,q]);
                    end
                end
            else
                fprintf('error CH_HL: v0 and v1 cannot both be 0.\n');
            end
        end

    % step 3: compute w_0, s_a, h_b, c and update w, Uc, Uh accordingly
        y_q =  bitget(y,q);
        %fprintf('y_q:');disp(dec2bin(y_q));
        s_a = uint8(0); 
        h_b = uint8(0);
        c = uint8(0);
        w_w = double(1);
        v_q = bitget(stab.v,q);
        %fprintf('v_q:');disp(dec2bin(v_q));
        if ~v_q
            if ~y_q
                if d == 0
                    %fprintf('~v_q, ~y_q, d==0\n');
                    s_a = 0; h_b = 1; c = 0; w_w = 2.^(0.5);
                elseif d == 1
                    s_a = 1; h_b = 1; c = 0; w_w = 2.^(0.5);
                elseif d == 2
                    s_a = 0; h_b = 1; c = 1; w_w = 2.^(0.5);
                elseif d == 3
                    s_a = 1; h_b = 1; c = 1; w_w = 2.^(0.5);
                else
                    fprintf('error CH_HL: invalid d.\n');
                end
            else
                if d == 0
                    s_a = 0; h_b = 1; c = 0; w_w = 2.^(0.5);
                elseif d == 1
                    s_a = 1; h_b = 1; c = 0; w_w = (-1i) * 2.^(0.5);
                elseif d == 2
                    s_a = 0; h_b = 1; c = 0; w_w = -(2.^(0.5));
                elseif d == 3
                    s_a = 1; h_b = 1; c = 0; w_w = (1i) * 2.^(0.5);
                else
                    fprintf('error CH_HL: invalid d.\n');
                end
            end
        else
            if ~y_q
                if d == 0
                    s_a = 0; h_b = 0; c = 0; w_w = 2.^(0.5);
                elseif d == 1
                    s_a = 1; h_b = 1; c = 1; w_w = 1 + 1i;
                elseif d == 2
                    s_a = 0; h_b = 0; c = 1; w_w = 2.^(0.5);
                elseif d == 3
                    s_a = 1; h_b = 1; c = 0; w_w = 1 - 1i;
                else
                    fprintf('error CH_HL: invalid d.\n');
                end
            else
                if d == 0
                    s_a = 0; h_b = 0; c = 0; w_w = 2.^(0.5);
                elseif d == 1
                    s_a = 1; h_b = 1; c = 0; w_w = 1 + 1i;
                elseif d == 2
                    s_a = 0; h_b = 0; c = 1; w_w = -(2.^(0.5));
                elseif d == 3
                    s_a = 1; h_b = 1; c = 1; w_w = 1 - 1i;
                else
                    fprintf('error CH_HL: invalid d.\n');
                end
            end
        end
        if s_a
            stab.CH_gate('SR',q);
        end
        stab.set_v(q,h_b);
        stab.s = y;
        stab.set_s(q,c);
        stab.w = stab.w * 2.^(-0.5) * ((-1).^double(a)) * w_w;
    else
        stab.w = stab.w * 2.^(-0.5) * ((-1).^double(a) + ((1i).^double(stab.get_g(p))) * ((-1).^double(b)));
        stab.s = t;
    end
    stab_new = stab;
end
