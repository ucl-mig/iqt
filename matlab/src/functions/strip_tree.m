% Remove some unnecessary data from the tree. 
function tree = strip_tree(tree)

    for idx = 1:length(tree)
         
        % You don't need linear transform for internal nodes:
        if not(tree{idx}.IsLeaf)
        tree{idx} = rmfield(tree{idx},'Transformation');
        end
        
        % You don't need 
        if tree{idx}.IsLeaf
        tree{idx} = rmfield(tree{idx},'SplitFeatureIndex');
        tree{idx} = rmfield(tree{idx},'FeatureThreshold');
        tree{idx} = rmfield(tree{idx},'LeftChildIndex');
        tree{idx} = rmfield(tree{idx},'RightChildIndex');
        end

        % Remove indices from all nodes.
        tree{idx} = rmfield(tree{idx},'TrainingDataIndices');
        tree{idx} = rmfield(tree{idx},'TestDataIndices');
        tree{idx} = rmfield(tree{idx},'InfoContent');
        tree{idx} = rmfield(tree{idx},'Visited');    
    end    
    
end