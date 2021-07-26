%% typedef without Simulink
classdef const
    properties (Constant)
        init_max_qubits = 8;
        typecast_str = append('uint',num2str(const.init_max_qubits));
        typecast_str_d = append('uint',num2str(const.init_max_qubits*2));
        init_uint = cast(0,const.typecast_str);
        init_duint = cast(0,const.typecast_str_d);
    end
end
