%%
%%
%% globals
%% CAVEAT!!!! Clear workspace before each test run!
len = 2;
vec_len = 2.^len;
num_gates = 1000;
rng('default');

%% init ch zero state
s = CH_state(len);
s.CH_init('zero');

gate_record = zeros(num_gates,3);

%% apply gates to construct s and test conjugate amplitude
for i = 1:num_gates
    fprintf('%dth test\n',i);
    gate_choice = randi(4,1,1);
    if gate_choice == 1  %SLs
        bit_choice = randi(len,1,1);
        s.CH_gate('SL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),-1];
    elseif gate_choice == 2  %HLs
        bit_choice = randi(len,1,1);
        s.CH_gate('HL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),-1];
    elseif gate_choice == 3 %CXL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        s.CH_gate('CXL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),bit_choice(2)];
    elseif gate_choice == 4 %CZL
        bit_choice = randperm(len,2); % use 1st as control and 2nd as result
        s.CH_gate('CZL',bit_choice);
        gate_record(i,:) = [gate_choice,bit_choice(1),bit_choice(2)];
    else
        fprintf('error invalid gate choice.\n');
    end
    s_conj = CH_state(len);
    s_conj.transpose(s);
    fprintf('print s:\n');
    s.pp_CH('ch');
    fprintf('print s_conj:\n');
    s.pp_CH('conj');
    assert(approx_equal(CH_CH_inner_product(s,s),1,0.000000001));
    %s_conj.pp_CH('ch');
    %s_conj.pp_CH('basis');
    fprintf('%dth test passed\n',i);
end

fprintf('conjugate test passed!\n');
