function decodedData = Huffman_decode(encodedData, root)

    decodedData = [];
    currentNode = root;

    for i = 1:length(encodedData)
        if encodedData(i) == '0'
            currentNode = currentNode.left;
        else
            currentNode = currentNode.right;
        end

        if isempty(currentNode.left) && isempty(currentNode.right)
            decodedData = [decodedData, currentNode.data];
            currentNode = root;
        end
    end
end