%%
%% class definition for stabilizer state in CH-form
%%
classdef CH_state < handle
    properties
        % caveat: all matrices are filled starting from top right (0,0) position

        % length/number of qubits of state
        len (1,1) int32 
        % len x len bit matrix storing Uc^*XpUc stabilizer tableau (X part)
        F = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % len x len bit matrix storing Uc^*ZpUc stabilizer tableau
        G = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % len x len bit matrix storing Uc^*XpUc stabilizer tableau (Z part)
        M = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % 2-bit x len array storing phase for Uc^*XpUc stabilizer tableau (phase part)
        g (1,1) = const.init_duint
        % len bit array storing Uh
        % len x len bit matrix storing UcXpUc^* stabilizer tableau (X part)
        FT = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % len x len bit matrix storing UcZpUc^* stabilizer tableau
        GT = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % len x len bit matrix storing UcXpUc^* stabilizer tableau (Z part)
        MT = cast(zeros(const.init_max_qubits,1),const.typecast_str)
        % 2-bit x len array storing phase for UcXpUc^* stabilizer tableau (phase part)
        gT (1,1) = const.init_duint
        % len bit array storing Uh
        v (1,1) = const.init_uint
        % len bit array storing |s>
        s (1,1) = const.init_uint
        % global phase on state
        w (1,1) double = 0
    end
%%------------------------------------------------------------------------------------------------------%%
    methods
        % initialize state given init option
        init(obj,option)
 
        % pretty print state given format option
        pp_CH(obj,option)
 
        % apply gate to state and return new state given gate g and control/target positions q
        stab_new = CH_gate(obj,g,q)
 
        % CH_state constructor: make meaningless stub state; don't forget to initialize with init
        function obj = CH_state(len) 
            obj.len = len;
        end

        % deepcopy: deepy copy the given state
        function deepcopy(obj,stab)
            obj.F = stab.F;
            obj.G = stab.G;
            obj.M = stab.M;
            obj.g = stab.g;
            obj.FT = stab.FT;
            obj.GT = stab.GT;
            obj.MT = stab.MT;
            obj.gT = stab.gT;
            obj.v = stab.v;
            obj.s = stab.s;
            obj.w = stab.w;
        end

        % transpose: make obj the transpose of stab without mutating stab
        function transpose(obj,stab)
            obj.F = stab.FT;
            obj.G = stab.GT;
            obj.M = stab.MT;
            obj.g = stab.gT;
            obj.FT = stab.F;
            obj.GT = stab.G;
            obj.MT = stab.M;
            obj.gT = stab.g;
            obj.v = stab.v;
            obj.s = stab.s;
            obj.w = conj(stab.w);
        end
 
        % CH_state property get and set
        % not required as properties are public, but this abstracts out bit ops and simplifies refactoring

        %getter
        function F_ij = get_F(obj,i,j) 
            F_ij = bitget(obj.F(i,1),j);
        end
       
        function G_ij = get_G(obj,i,j) 
            G_ij = bitget(obj.G(i,1),j);
        end
 
        function M_ij = get_M(obj,i,j) 
            M_ij = bitget(obj.M(i,1),j);
        end
 
        function g_i = get_g(obj,i) 
            g_i1 = bitget(obj.g,2*i-1);
            g_i2 = bitshift(bitget(obj.g,2*i),1);
            g_i = bitxor(g_i1,g_i2);
            g_i = cast(g_i,'like',const.init_uint);

        end

        function FT_ij = get_FT(obj,i,j) 
            FT_ij = bitget(obj.FT(i,1),j);
        end
       
        function GT_ij = get_GT(obj,i,j) 
            GT_ij = bitget(obj.GT(i,1),j);
        end
 
        function MT_ij = get_MT(obj,i,j) 
            MT_ij = bitget(obj.MT(i,1),j);
        end
 
        function gT_i = get_gT(obj,i) 
            gT_i1 = bitget(obj.gT,2*i-1);
            gT_i2 = bitshift(bitget(obj.gT,2*i),1);
            gT_i = bitxor(gT_i1,gT_i2);
            gT_i = cast(gT_i,'like',const.init_uint);
        end
 
        function v_i = get_v(obj,i) 
            v_i = bitget(obj.v,i);
        end
 
        function s_i = get_s(obj,i) 
            s_i = bitget(obj.s,i);
        end
       
        % setter
        function set_F(obj,i,j,b) 
            obj.F(i,1) = bitset(obj.F(i,1),j,b);
        end
 
        function set_G(obj,i,j,b) 
            obj.G(i,1) = bitset(obj.G(i,1),j,b);
        end
 
        function set_M(obj,i,j,b) 
            obj.M(i,1) = bitset(obj.M(i,1),j,b);
        end
 
        function set_g(obj,i,b) 
            obj.g = bitset(obj.g,2*i-1,bitget(b,1));
            obj.g = bitset(obj.g,2*i,bitget(b,2));
        end

        function set_FT(obj,i,j,b) 
            obj.FT(i,1) = bitset(obj.FT(i,1),j,b);
        end
 
        function set_GT(obj,i,j,b) 
            obj.GT(i,1) = bitset(obj.GT(i,1),j,b);
        end
 
        function set_MT(obj,i,j,b) 
            obj.MT(i,1) = bitset(obj.MT(i,1),j,b);
        end
 
        function set_gT(obj,i,b) 
            obj.gT = bitset(obj.gT,2*i-1,bitget(b,1));
            obj.gT = bitset(obj.gT,2*i,bitget(b,2));
        end
 
        function set_v(obj,i,b) 
            obj.v = bitset(obj.v,i,b);
        end

        function set_s(obj,i,b) 
            obj.s = bitset(obj.s,i,b);
        end
    end
end
