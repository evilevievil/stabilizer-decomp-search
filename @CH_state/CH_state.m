%%
%% class definition for len-qubit stabilizer state in CH-form
%%
classdef CH_state
   properties
      % todo: max 8 qubits for now. make global param/improve readability later
      % all matrices are filled starting from top right (0,0) position
      % length/number of qubits of state
      len (1,1) int32 {mustBeNonnegative,mustBeInteger, mustBeFinite} 
      % len x len bit matrix storing Uc^*ZpUc stabilizer tableau
      F (8,1) uint8 = zeros(8,1)
      % len x len bit matrix storing Uc^*XpUc stabilizer tableau (X part)
      G (8,1) uint8 = zeros(8,1)
      % len x len bit matrix storing Uc^*XpUc stabilizer tableau (Z part)
      M (8,1) uint8 = zeros(8,1)
      % 4-bit x len array storing phase for Uc^*XpUc stabilizer tableau (phase part)
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

      % computes inner product <a|obj>
      inner_product = CH_inner_product(obj,a)

      % CH_state constructor: make meaningless stub state; don't forget to initialize with init
      function obj = CH_state(len) 
         obj.len = len;
      end

      % CH_state property get and set
      % not required as properties are public, but this abstracts out bit ops and simplifies refactoring
      % getter
      function F_ij = get.F(obj,i,j) 
         F_ij = bitget((obj.F)(i,1),j);
      end

      function G_ij = get.G(obj,i,j) 
         G_ij = bitget((obj.G)(i,1),j);
      end

      function M_ij = get.M(obj,i,j) 
         M_ij = bitget((obj.M)(i,1),j);
      end

      function g_i = get.g(obj,i) 
         g_i1 = uint8(bitget(obj.g,2*i-1));
         g_i2 = bitshift(uint8(bitget(obj.g,2*i)),1);
         g_i = bitand(g_i1,g_i2);
      end

      function v_i = get.v(obj,i) 
         v_i = bitget(obj.v,i);
      end

      function s_i = get.s(obj,i) 
         s_i = bitget(obj.s,i);
      end

      function w = get.w(obj) 
         w = obj.w;
      end
      
      % setter
      function set.F(obj,i,j,b) 
         bitset((obj.F)(i,1),j,b);
      end

      function set.G(obj,i,j,b) 
         bitset((obj.G)(i,1),j,b);
      end

      function set.M(obj,i,j,b) 
         bitset((obj.M)(i,1),j,b);
      end

      function set.g(obj,i,b) 
         bitset(obj.g,2*i-1,bitget(b,1));
         bitset(obj.g,2*i,bitget(b,2));
      end

      function set.v(obj,i,b) 
         bitset(obj.v,i,b);
      end

      function set.s(obj,i,b) 
         bitset(obj.s,i,b);
      end

      function set.w(obj,val) 
         obj.w = val;
      end
   end
end