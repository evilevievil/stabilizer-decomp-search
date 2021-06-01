%%
%% class definition for len-qubit stabilizer state in CH-form
%%
classdef CH_state < handle
    properties
        % todo: max 8 qubits for now. make global param/improve readability later
        % todo: make functions more robust by asserting unuused bit equal 0
        % all matrices are filled starting from top right (0,0) position
        % length/number of qubits of state
        len (1,1) int32 {mustBeNonnegative,mustBeInteger, mustBeFinite} 
        % len x len bit matrix storing Uc^*XpUc stabilizer tableau (X part)
        F (8,1) uint8 = zeros(8,1)
        % len x len bit matrix storing Uc^*ZpUc stabilizer tableau
        G (8,1) uint8 = zeros(8,1)
        % len x len bit matrix storing Uc^*XpUc stabilizer tableau (Z part)
        M (8,1) uint8 = zeros(8,1)
        % 2-bit x len array storing phase for Uc^*XpUc stabilizer tableau (phase part)
        g (1,1) uint16 = 0
        % len bit array storing Uh
        v (1,1) uint8 = 0
        % len bit array storing |s>
        s (1,1) uint8 = 0
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
            obj.v = stab.v;
            obj.s = stab.s;
            obj.w = stab.w;
        end

        % transpose: make obj the transpose of stab without mutating stab
        function transpose(obj,stab)
            %obj.F = stab.F;
            %obj.G = stab.G;
            %obj.M = stab.M;
            %obj.g = stab.g;
            obj.v = stab.v;
            obj.s = stab.s;
            obj.w = stab.w;
        end
 
        % CH_state property get and set
        % not required as properties are public, but this abstracts out bit ops and simplifies refactoring
        % getter
        % todo: make inline/anaonymous functions to boost performance 
 
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
            g_i1 = uint8(bitget(obj.g,2*i-1));
            g_i2 = bitshift(uint8(bitget(obj.g,2*i)),1);
            g_i = bitxor(g_i1,g_i2);
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
 
        function set_v(obj,i,b) 
            obj.v = bitset(obj.v,i,b);
        end

        function set_s(obj,i,b) 
            obj.s = bitset(obj.s,i,b);
        end
    end
end