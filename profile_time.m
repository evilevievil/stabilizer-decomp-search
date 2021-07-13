%b = magic_state_vec('catT',3);
%a = kron(b,b);
%fixed_rank_stab_decomp_search(a,6,4,200,200,1000,1);

%b = magic_state_vec('catT',3);
%a = kron(b,b);
%fixed_rank_stab_decomp_search(b,3,4,200,200,1000,1);

%a = magic_state_vec('H',6);
%fixed_rank_stab_decomp_search(a,6,6,400,4000,100,1000);
%fixed_rank_stab_decomp_search(a,5,5,1,300,100,1000); % H_5_5_0.9958
%fixed_rank_stab_decomp_search(a,6,7,1,300,100,1000); % H_6_7_0.9539 pt1
%H_6_7_0.9539 pt2 -> 0.9898 4000
%% current
a = magic_state_vec('T',6);
%fixed_rank_stab_decomp_search(a,6,7,1,4000,100,1000); % H_6_7_0.9539 pt2 -> 0.9898
fixed_rank_stab_decomp_search(a,6,6,1,4000,100,1000); % H_6_7_0.9539 pt2 -> 0.9898

% a = magic_state_vec('H',6);
% rng(89);  
% for i = 1:100
%   seed = randi(1000000000,1,1);
%   fixed_rank_stab_decomp_search(a,6,7,1,4000,200,100,seed);
% end










