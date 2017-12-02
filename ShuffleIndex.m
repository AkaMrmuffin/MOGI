function [ c_sorted, indices] = ShuffleIndex( c_sorted, indices )

L = length(c_sorted);
indices_needing_shuffling = [];
for i = 2:1:L
    if c_sorted(i) == c_sorted (i-1)
        catindex = cat(1,indices_needing_shuffling,i);
        indices_needing_shuffling = [indices_needing_shuffling catindex];    
   
        if length(indices_needing_shuffling)>0
           shuffled_indeces = randsample(indices_needing_shuffling, length(indices_needing_shuffling));

           for j = 1:length(indices_needing_shuffling)
               new_i = indices_needing_shuffling(j);
               indices (new_i) = shuffled_indeces(j);
           end
           indices_needing_shuffling = [];
        end
    end
end


end

