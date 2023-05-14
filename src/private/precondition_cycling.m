function [] = precondition_cycling(array, index, strategy, memory)

[isrv, ~]  = isrealvector(array);
assert(isrv);  

assert(isintegerscalar(index));

assert(isintegerscalar(strategy) && 0<=strategy && strategy<=5);

assert(islogicalscalar(memory));
end
