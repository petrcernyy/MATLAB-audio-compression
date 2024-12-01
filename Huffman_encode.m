function [root, encoded_data] = Huffman_encode(data)

    unique_data = unique(data);
    character_freq = zeros(1,length(unique_data));
    for i = 1:length(unique_data)
        character_freq(i) = length(find(unique_data(i)==data));
    end

    n = length(unique_data);
    nodes(n) = newNode([], []);
    for i = 1:n
        nodes(i) = newNode(unique_data(i), character_freq(i));
    end
    
    [~, idx] = sort([nodes.freq], 'ascend');
    nodes = nodes(idx);

    while numel(nodes) > 1
        left = nodes(1);
        nodes(1) = [];

        right = nodes(1);
        nodes(1) = [];

        top = newNode('-', left.freq + right.freq);
        top.left = left;
        top.right = right;

        nodes(end+1) = top;

        [~, idx] = sort([nodes.freq], 'ascend');
        nodes = nodes(idx);
    end

    root = nodes(1);

    global counter;
    global new_freq;
    
    new_freq = cell(length(character_freq),2);
    counter = 0;

    iterative_print_codes(root, '');

    encoded_data = '';
    for i = 1:length(data)
        encoded_data = append(encoded_data, new_freq{find([new_freq{:,1}] == data(i)),2});
    end

end

function node = newNode(data, freq)
    node.data = data;
    node.freq = freq;
    node.left = [];
    node.right = [];
end

function iterative_print_codes(root, code)

    global counter;
    global new_freq;

    if isempty(root.left) && isempty(root.right)
        counter = counter + 1;
        new_freq{counter,1} = root.data;
        new_freq{counter,2} = code;
    else
        iterative_print_codes(root.left, strcat(code, '0'));
        iterative_print_codes(root.right, strcat(code, '1'));
    end

end
