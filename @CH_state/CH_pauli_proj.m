%%
%%  stab_new = stab.CH_pauli_proj(g,q)
%%               
function stab_new = CH_pauli_proj(obj,g,q)
    switch g
    case '+Z'
        stab_new = CH_Z(0,q,obj);
    case '-Z'
        stab_new = CH_Z(1,q,obj);
    case '+X'
        stab_new = CH_X(0,q,obj);
    case '-X'
        stab_new = CH_X(1,q,obj);
    case '+Y'
        stab_new = CH_Y(0,q,obj);
    case '-Y'
        stab_new = CH_Y(1,q,obj);
    otherwise
        fprintf('error CH_pauli_proj: invalid pauli projector name.\n');
    end
end


function stab_new = CH_Z(neg,q,stab)
    stab_new = CH_state(stab.len);
    stab_new.deepcopy(stab);
    Zq = stab_new.G(q,1);
    t = stab.s;
    u = bitxor(bitand(Zq,stab_new.v),t); % gives X positions
    a = uint8(0);
    d = uint8(0); %no additional phase for Z
    neg = parity(neg+bitand(bitxor(bitand(Zq,stab_new.v),Zq),t));
    if neg
        d = 2;
    end

    if t ~= u
        % todo: refactor superpos2circuit to take 2.^(0.5) out 
        superpos2circuit(t,u,a,d,stab_new);
        stab_new.w = 0.5 * stab_new.w * 2.^(0.5);
    else
        if neg
            stab_new.w = 0;
        end
    end
end

function stab_new = CH_X(neg,q,stab)
    stab_new = CH_state(stab.len);
    stab_new.deepcopy(stab);
    %second layer
    ZqL = bitand(stab_new.F(q,1),stab_new.v);
    XqR = bitand(stab_new.M(q,1),stab_new.v);
    % first layer
    ZqR = bitxor(XqR,stab_new.M(q,1));
    XqL = bitxor(ZqL,stab_new.F(q,1));
    t = stab.s;
    neg = parity(neg+bitand(ZqR,t));
    u = bitxor(XqR,t); 
    neg = parity(neg+bitand(ZqL,u));
    u = bitxor(XqL,u); 
    a = uint8(0);
    d = stab_new.get_g(q);
    if neg
        d = bitset(d,2,~bitget(d,2));
    end

    if t ~= u
        % todo: refactor superpos2circuit to take 2.^(0.5) out 
        superpos2circuit(t,u,a,d,stab_new);
        stab_new.w = 0.5 * stab_new.w * 2.^(0.5);
    else
        stab_new.w = stab_new.w * 0.5 * (1+(1i).^double(d));
    end
end

function stab_new = CH_Y(neg,q,stab)
    stab_new = CH_state(stab.len);
    stab_new.deepcopy(stab);
    % zeroth layer
    Zq = stab_new.G(q,1);
    Xq = bitand(Zq,stab_new.v);
    t = stab.s;
    u = bitxor(Xq,t); % gives X positions
    a = uint8(0);
    d = mod(stab_new.get_g(q)+1,4);
    neg = parity(neg+bitand(bitxor(bitand(Zq,stab_new.v),Zq),t));
    %second layer
    ZqL = bitand(stab_new.F(q,1),stab_new.v);
    XqR = bitand(stab_new.M(q,1),stab_new.v);
    % first layer
    ZqR = bitxor(XqR,stab_new.M(q,1));
    XqL = bitxor(ZqL,stab_new.F(q,1));
    neg = parity(neg+bitand(ZqR,u));
    u = bitxor(XqR,u); 
    neg = parity(neg+bitand(ZqL,u));
    u = bitxor(XqL,u); 
    if neg
        d = bitset(d,2,~bitget(d,2));
    end

    if t ~= u
        % todo: refactor superpos2circuit to take 2.^(0.5) out 
        fprintf('superpos case\n');
        superpos2circuit(t,u,a,d,stab_new);
        stab_new.w = 0.5 * stab_new.w * 2.^(0.5);
    else
        stab_new.w = stab_new.w * 0.5 * (1+(1i).^double(d));
    end
end
