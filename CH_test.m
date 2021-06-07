%%
%% tests/examples for CH_state methods/functions
%%

%% initialization sanity checks
s = CH_state(4);
s.CH_init('zero');
%s.pp_CH('ch');
assert(isequal(s.F,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.G,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.M,[0;0;0;0;0;0;0;0]));
assert(s.g == uint16(0));
assert(s.v == uint8(0));
assert(s.s == uint8(0));
assert(s.w == double(1));
fprintf('init sanity test passed!\n');

%% CH_gate sanity checks

% SL
s = CH_state(4);
s.CH_init('zero');
s.CH_gate('SL',3);
assert(isequal(s.M,[0;0;bin2dec('0100');0;0;0;0;0]));
assert(s.g == uint16(bin2dec('110000')));
s.CH_gate('SL',3);
assert(isequal(s.M,[0;0;0;0;0;0;0;0]));
assert(s.g == uint16(bin2dec('100000')));
s.CH_gate('SL',3);
assert(isequal(s.M,[0;0;bin2dec('0100');0;0;0;0;0]));
assert(s.g == uint16(bin2dec('10000')));
s.CH_gate('SL',3);
assert(isequal(s.M,[0;0;0;0;0;0;0;0]));
assert(s.g == uint16(bin2dec('0')));
fprintf('SL gate sanity test passed!\n');

% SR
s = CH_state(4);
s.CH_init('zero');
s.CH_gate('SR',2);
assert(isequal(s.M,[0;bin2dec('0010');0;0;0;0;0;0]));
assert(s.g == uint16(bin2dec('1100')));
s.set_F(2,2,0); 
s.CH_gate('SR',2);
assert(isequal(s.M,[0;bin2dec('0010');0;0;0;0;0;0]));
assert(s.g == uint16(bin2dec('1100')));
fprintf('SR gate sanity test passed!\n');

% CZR
s = CH_state(4);
s.M(1,1) = bin2dec('1010');
s.M(2,1) = bin2dec('0010');
s.M(3,1) = bin2dec('1001');
s.M(4,1) = bin2dec('0011');
s.F(1,1) = bin2dec('1011');
s.F(2,1) = bin2dec('0111');
s.F(3,1) = bin2dec('1011');
s.F(4,1) = bin2dec('0011');
s.g = uint16(bin2dec('11100100'));
s.CH_gate('CZR',[1,2]);
assert(isequal(s.M,[bin2dec('1001'); ...
                    bin2dec('0001'); ...
                    bin2dec('1010'); ...
                    bin2dec('0000'); ...
                    0;0;0;0]));
assert(s.g == uint16(bin2dec('01001110')));
s.CH_gate('CZR',[3,4]);
assert(isequal(s.M,[bin2dec('1101'); ...
                    bin2dec('1001'); ...
                    bin2dec('1110'); ...
                    bin2dec('0000'); ...
                    0;0;0;0]));
assert(s.g == uint16(bin2dec('01001110')));
fprintf('CZR gate sanity test passed!\n');

% CZL
s.G(1,1) = bin2dec('1010');
s.G(2,1) = bin2dec('0010');
s.G(3,1) = bin2dec('1001');
s.G(4,1) = bin2dec('0011');
s.CH_gate('CZL',[3,4]);
assert(isequal(s.M,[bin2dec('1101'); ...
                    bin2dec('1001'); ...
                    bin2dec('1101'); ...
                    bin2dec('1001'); ...
                    0;0;0;0]));
fprintf('CZL gate sanity test passed!\n');

% CXL
s = CH_state(4);
s.M(1,1) = bin2dec('1010');
s.M(2,1) = bin2dec('0010');
s.M(3,1) = bin2dec('1001');
s.M(4,1) = bin2dec('0011');
s.F(1,1) = bin2dec('1011');
s.F(2,1) = bin2dec('0111');
s.F(3,1) = bin2dec('1011');
s.F(4,1) = bin2dec('0011');
s.G(1,1) = bin2dec('1010');
s.G(2,1) = bin2dec('0010');
s.G(3,1) = bin2dec('1001');
s.G(4,1) = bin2dec('0011');
s.CH_gate('CXL',[3,2]);
s.CH_gate('CXL',[4,1]);
assert(isequal(s.F,[bin2dec('1011'); ...
                    bin2dec('0111'); ...
                    bin2dec('1100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.G,[bin2dec('1001'); ...
                    bin2dec('1011'); ...
                    bin2dec('1001'); ...
                    bin2dec('0011'); ...
                    0;0;0;0]));
assert(isequal(s.M,[bin2dec('1010'); ...
                    bin2dec('0010'); ...
                    bin2dec('1011'); ...
                    bin2dec('1001'); ...
                    0;0;0;0])); 
s.g = uint16(bin2dec('00100100'));
s.CH_gate('CXL',[3,2]);
s.CH_gate('CXL',[4,1]);
assert(s.g == uint16(bin2dec('00110100')));
fprintf('CXL gate sanity test passed!\n');

% CXR
s = CH_state(4);
s.M(1,1) = bin2dec('1010');
s.M(2,1) = bin2dec('0010');
s.M(3,1) = bin2dec('1001');
s.M(4,1) = bin2dec('0011');
s.F(1,1) = bin2dec('1011');
s.F(2,1) = bin2dec('0111');
s.F(3,1) = bin2dec('1011');
s.F(4,1) = bin2dec('0011');
s.G(1,1) = bin2dec('1010');
s.G(2,1) = bin2dec('0010');
s.G(3,1) = bin2dec('1001');
s.G(4,1) = bin2dec('0011');
s.CH_gate('CXR',[3,2]);
s.CH_gate('CXR',[4,1]);
assert(isequal(s.F,[bin2dec('1010'); ...
                    bin2dec('0101'); ...
                    bin2dec('1010'); ...
                    bin2dec('0011'); ...
                    0;0;0;0]));
assert(isequal(s.G,[bin2dec('1110'); ...
                    bin2dec('0110'); ...
                    bin2dec('0001'); ...
                    bin2dec('1111'); ...
                    0;0;0;0]));
assert(isequal(s.M,[bin2dec('1110'); ...
                    bin2dec('0110'); ...
                    bin2dec('0001'); ...
                    bin2dec('1111'); ...
                    0;0;0;0]));                
fprintf('CXR gate sanity test passed!\n');

% HL
% simple HL test 1
s = CH_state(4);
s.CH_init('zero');
s.CH_gate('HL',1);
s.CH_gate('HL',2);
s.CH_gate('HL',4);
assert(isequal(s.F,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.G,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.M,[0;0;0;0;0;0;0;0]));
assert(s.g == uint16(0));
assert(s.v==bin2dec('1011'));
assert(s.s == uint8(0));
assert(s.w == double(1));

% simple HL test 2
s = CH_state(4);
s.CH_init('zero');
s.CH_gate('HL',1);
s.CH_gate('SL',1);
s.CH_gate('HL',1);
assert(isequal(s.F,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.G,[bin2dec('0001'); ...
                    bin2dec('0010'); ...
                    bin2dec('0100'); ...
                    bin2dec('1000'); ...
                    0;0;0;0]));
assert(isequal(s.M,[bin2dec('0001');0;0;0;0;0;0;0]));
assert(s.g == uint16(bin2dec('11')));
assert(s.v==bin2dec('0001'));
assert(s.s == bin2dec('0001'));
assert(s.w == (1+1i) * 2.^(-0.5));
fprintf('HL gate simple sanity tests passed!\n');

% todo: sanity tests for pauli projector and conjugate?

%% inner product/projection sanity tests

% basis stab inner product sanity check
s = CH_state(2);
s.CH_init('zero');
s.CH_gate('HL',1);
assert(CH_basis_inner_product(bin2dec('00'),s)==2.^(-0.5));
assert(CH_basis_inner_product(bin2dec('01'),s)==2.^(-0.5));
assert(CH_basis_inner_product(bin2dec('10'),s)==0);
assert(CH_basis_inner_product(bin2dec('11'),s)==0);
s.CH_gate('SL',1);
assert(CH_basis_inner_product(bin2dec('00'),s)==2.^(-0.5));
assert(CH_basis_inner_product(bin2dec('01'),s)==2.^(-0.5)*(1i));
assert(CH_basis_inner_product(bin2dec('10'),s)==0);
assert(CH_basis_inner_product(bin2dec('11'),s)==0);
s.CH_gate('CXL',[1,2]);
assert(CH_basis_inner_product(bin2dec('00'),s)==2.^(-0.5));
assert(CH_basis_inner_product(bin2dec('01'),s)==0);
assert(CH_basis_inner_product(bin2dec('10'),s)==0);
assert(CH_basis_inner_product(bin2dec('11'),s)==2.^(-0.5)*(1i));
fprintf('basis stab inner product tests passed!\n');

% stab stab inner product sanity check
s1 = CH_state(2);
s1.CH_init('zero');
s1.CH_gate('HL',1);

s2 = CH_state(2);
s2.CH_init('zero');
s2.CH_gate('HL',1);

assert(CH_CH_inner_product(s1,s2)==1);
s2.CH_gate('HL',2);
assert(CH_CH_inner_product(s1,s2)==2.^-0.5);
s2.CH_gate('SL',1);
assert(approx_equal(CH_CH_inner_product(s1,s2),(1+1i)*(2.^-1.5),0.000000001));
assert(approx_equal(CH_CH_inner_product(s2,s1),(1-1i)*(2.^-1.5),0.000000001));
fprintf('stab stab inner product tests passed!\n');

% Projection Sanity tests
s3 = CH_state(2);
s3.CH_init('zero');
s3.CH_gate('HL',1);
a = zeros(1,4); % (1/root2,0,1/root2,0)
a(1) = 2.^-0.5;
a(2) = 2.^-0.5; %follow reverse string convention 
a(3) = sqrt(0);  %follow reverse string convention 
a(4) = sqrt(0);
disp(CH_decomp_project(a,[s3],2,1));
assert(approx_equal(CH_decomp_project(a,[s3],2,1),1,0.000000001));
fprintf('stab decomp projection tests passed!\n');
