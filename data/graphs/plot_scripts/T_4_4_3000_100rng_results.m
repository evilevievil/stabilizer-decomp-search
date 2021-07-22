
num_sucess=zeros(1,25);
total_num_sucess = 0;
a = magic_state_vec('T',3);
for i =1:25
   rng(i);
   [success_status,decomp] = fixed_rank_stab_decomp_search(a,3,3,1,6000,100,1000);
   num_sucess(1,i) = success_status;
   total_num_sucess = total_num_sucess + success_status;
end

disp(total_num_sucess);
