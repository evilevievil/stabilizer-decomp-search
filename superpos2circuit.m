%% 
%%  superpos2circuit(t,u,a,d,stab)
%%  absorb superpoisition state (assumes t!=u) into circuit and mutates stab
%%  (see prop'n 4 from "Simulation of qc by low-rank stab decomp" paper)
%% 
function superpos2circuit(t,u,a,d,stab)
    % compute strings y,z (differ by one bit q), string v0,v1 representing Vc
    t_oxr_u = bitxor(t,u);
    v0 = bitand(t_oxr_u,bitcmp(stab.v));
    v1 = bitand(t_oxr_u,stab.v);
    q = 1;

    if v0
        while ~bitget(v0,q)
            q = q+1;
        end
    elseif v1
        while ~bitget(v1,q)
            q = q+1;
        end
    else
        fprintf('error CH_HL: v0 and v1 cannot both be 0.\n');
    end
    v_q = bitget(stab.v,q);

    % update Uc = UcVc
    if bitsum(t_oxr_u) > 1
        if v0 
            for pos = 1:stab.len
                if bitget(v0,pos) && (pos ~= q)
                    stab.CH_gate('CXR',[q,pos]);
                end
            end
            for pos = 1:stab.len
                if bitget(v1,pos)
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

    if bitget(t,q)
        y = bitxor(u,bitset(const.init_uint,q));
    else
        y = t;
    end

    % compute w_0, s_a, h_b, c and update w, Uc, Uh accordingly
    y_q =  bitget(y,q);
    s_a = const.init_uint; 
    h_b = const.init_uint;
    c = const.init_uint;
    w_w = double(1);
    
    % todo: can simplify case analysis
    if ~v_q
        if ~y_q
            if d == 0
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
                s_a = 1; h_b = 1; c = 1; w_w = (1i) * 2.^(0.5);
            elseif d == 2
                s_a = 0; h_b = 1; c = 1; w_w = -(2.^(0.5));
            elseif d == 3
                s_a = 1; h_b = 1; c = 0; w_w = (-1i) * 2.^(0.5);
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
end
