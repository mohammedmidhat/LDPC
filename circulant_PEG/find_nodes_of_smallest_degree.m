function node = find_nodes_of_smallest_degree(array, chk_degree_sequence)
    [sorted_array Indices] = sort(array,'ascend');
    for i=sorted_array
        if(array(1,i))
            node = i;
            break;
    end
end
