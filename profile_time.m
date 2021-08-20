%% random/sporadic tests and experiments 

% fixed-rank search for various magic states

%b = magic_state_vec('catT',3);
%a = kron(b,b);
%fixed_rank_stab_decomp_search(a,6,4,200,200,1000,1);

%a = magic_state_vec('H',6);
%fixed_rank_stab_decomp_search(a,6,6,400,4000,100,1000);

a = magic_state_vec('H',5);
fixed_rank_stab_decomp_search(a,5,5,1,300,100,1000); % H_5_5_0.9958

%a = magic_state_vec('T',6);
%fixed_rank_stab_decomp_search(a,6,7,1,4000,100,1000);

%a = magic_state_vec('r_1_3',8);
%fixed_rank_stab_decomp_search(a,8,2,200,3000,100,2000);

%a = magic_state_vec('12gencat',12);
%fixed_rank_stab_decomp_search(a,12,5,5000,5000,100,1000);


% magic code state search
%magic_state_decomp_search_v2(6,3,1,4000,200,1000,5,2,4); 
%magic_state_decomp_search_v2(8,3,1,4000,200,1000,5,-1,4); %0.2436
%magic_state_decomp_search_v2(8,3,1,4000,200,1000,5,2,4); %0.2436
