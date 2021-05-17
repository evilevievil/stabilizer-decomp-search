%%
%%  stab_new = CH_gate(g,q)
%%                                             
function stab_new = CH_gate(g,q)
end

% apply C-type gates from left
function stab_new = CH_SL(p,stab)
end

function stab_new = CH_CZL(p,stab)
end

function stab_new = CH_CXL(p,stab)
end

% apply C-type gates from right 
function stab_new = CH_SR(p,stab)
end

function stab_new = CH_CZR(p,stab)
end

function stab_new = CH_CXR(p,stab)
end

% apply H gate on the left
function stab_new = CH_H(p,stab)
end
