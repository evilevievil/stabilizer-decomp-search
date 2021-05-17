%% stab_decomp_search(a)
%% main algorithm that searches for low-rank stabilizer decomposition of arbitrary state a

function stab_decomp = stab_decomp_search(a)
    %%% init
    % rank k
    % upper limit on k RANK_UPPER_LIMIT
    % start state s (try 3 options: 1. user input; 2. rand start state 3. optimized start state)
    % ? check if a is a stabilizer state ?

    %%% main (SA algo)
    % for k=1 to RANK_UPPER_LIMIT do:
    %     upper limit on SA iterations SA_UPPER_LIMIT
    %     run SA until stab_decomp found or SA_UPPER_LIMIT is reached
    
    %% print summary
    %  summarize search results for each k in chart/graph
    %  store search results in db? could be used as training set?
    
end
