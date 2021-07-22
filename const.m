%% typedef without Simulink
classdef const
    properties (Constant)
        typecast_str = 'uint16';
        init_uint uint16 = 0;
        init_duint uint32 = 0;
        init_tableau (16,1) uint16 = zeros(const.init_max_qubits,1);
        init_max_qubits = 16;
    end
end
