function order = conn_statsoptimalleaforder(Z,Y)
%OPTIMALLEAFORDER optimal leaf ordering for hierarchical clustering.

%   Reference:
%   Bar-Joseph Z, Gifford DK, Jaakkola TS. Fast optimal leaf ordering for
%   hierarchical clustering. Bioinformatics. 2001;17 Suppl 1:S22-9. PMID:
%   11472989

% Copyright 2006-2012 The MathWorks, Inc.


numLeaves = size(Z,1)+1;
% check the tree
if size(Z,2) == 3
    Z(:,3)=[];
end
Y = (Y+Y')/2;
if numLeaves == 2
    order = [1 2];
    return;
end

groupCriteria = false;
transFunction = @(d) max(d(:))-d;

% transform to similarities
mask = true(1,numLeaves*numLeaves);
mask(1:numLeaves+1:numLeaves*numLeaves)=false;
Y(mask) = transFunction(Y(mask));

% initialize the scoring and traceback matrices
scores = -inf(numLeaves);
tracebackMatrix = zeros(numLeaves,'uint32');

% main loop
for outer = 1:size(Z,1)-1
    % get the leave of the node
    leftLeaves = getChildrenFromNode(Z(outer,1));
    rightLeaves = getChildrenFromNode(Z(outer,2));
    numLeft = numel(leftLeaves);
    numRight = numel(rightLeaves);
    dm = Y(rightLeaves,leftLeaves);
    if groupCriteria
       srdm=sum(dm,2); scdm=sum(dm,1); 
       dm = srdm(:,ones(numLeft,1))+scdm(ones(numRight,1),:);
    end
    if numLeft == 1
        if numRight == 1
            % If we are at a node that has just two leaves then there is
            % nothing to calculate and we only populate the traceback matrices.
            scores(leftLeaves,rightLeaves) = dm;
            scores(rightLeaves,leftLeaves) = dm;
            tracebackMatrix(rightLeaves,leftLeaves) = rightLeaves;
            tracebackMatrix(leftLeaves,rightLeaves) = leftLeaves;
        else
            % We are at a node that has one leaf on the left and a branch
            % on the right; use a vectorized approach to find a right-inner
            % leaf that leads to a maximum score for every possible
            % right-outer leaf:  
           % [BEST,TR] = max(dm(:,ones(numRight,1))+scores(rightLeaves,rightLeaves)); 
            [BEST,TR] = max(bsxfun(@plus,dm,scores(rightLeaves,rightLeaves))); 
            scores(leftLeaves,rightLeaves) = BEST;
            scores(rightLeaves,leftLeaves) = BEST;
            scores(rightLeaves,rightLeaves) = -inf; 
            %save right_inner leaf
            tracebackMatrix(rightLeaves,leftLeaves)  = rightLeaves(TR);
            tracebackMatrix(leftLeaves, rightLeaves) = leftLeaves;
        end
    elseif numRight == 1
        % We are at a node that has one leaf on the right and a branch on
        % the left; use a vectorized approach to find a left-inner leaf
        % that leads to a maximum score for every possible left-outer leaf: 
        [BEST,TL] = max(dm(ones(numLeft,1),:)'+scores(leftLeaves,leftLeaves));  
        scores(leftLeaves,rightLeaves) = BEST;
        scores(rightLeaves,leftLeaves) = BEST;
        scores(leftLeaves,leftLeaves) = -inf;
        tracebackMatrix(rightLeaves,leftLeaves)  = rightLeaves;
        tracebackMatrix(leftLeaves, rightLeaves) = leftLeaves(TL);
    else
        % We are at a node that joins two branches, each one with more than
        % one leave, sort the scores as we will loop over the highest
        % scores until we know that the next combination cannot be better
        % than the existing best (Bar-Joseph, et.al., 2001):
        [maxL,lPerm] = sort((scores(leftLeaves,leftLeaves)),1,'descend');
        [maxR,rPerm] = sort((scores(rightLeaves,rightLeaves)),1,'descend');
        
        dm_maxL = max(dm);
        dm_maxLR = max(dm_maxL);
        for left = 1:numLeft 
            for right = 1:numRight 
                
                best = -inf; 
                tL = 0;
                tR = 0;
                % set threshold to end global search
                maxPosR = maxR(1,right) + dm_maxLR;
                
                for iiL = 1:numLeft-1
                    innerL = lPerm(iiL,left); 
                    scoreL = maxL(iiL,left); 
                    if scoreL + maxPosR <= best
                        break
                    end
                    % set threshold to end search for best right-inner
                    % leaf, Bar-Joseph et.al., used dm_maxLR, but
                    % dm_maxL(innerL) further improves the speed up
                    maxPosL = scoreL + dm_maxL(innerL);
                    for iiR = 1:numRight-1
                        innerR = rPerm(iiR,right); 
                        scoreR = maxR(iiR,right); 
                        if maxPosL + scoreR <= best
                            break
                        end
                    
                        val = scoreL + scoreR + dm(innerR,innerL);
                        if val>best
                            best = val;
                            tL = leftLeaves(innerL);
                            tR = rightLeaves(innerR);
                        end
                    end % for iiR = 1:numRight-1
                end % for iiL = 1:numLeft-1
                % populate the traceback information
                scores(leftLeaves(left),rightLeaves(right)) = best;
                scores(rightLeaves(right),leftLeaves(left)) = best;
                tracebackMatrix(rightLeaves(right),leftLeaves(left)) = tR;
                tracebackMatrix(leftLeaves(left), rightLeaves(right)) = tL;
            end % for right = 1:numRight
        end % for left = 1:numLeft
        scores(leftLeaves,leftLeaves) = -inf;
        scores(rightLeaves,rightLeaves) = -inf;
    end % numLeft > 1  && numRight > 1
end % for outer = 1:size(Z,1)-1

%process the last branch separately to improve the efficiency
leftLeaves = getChildrenFromNode(Z(end,1));
rightLeaves = getChildrenFromNode(Z(end,2));
numLeft = numel(leftLeaves);
numRight = numel(rightLeaves);
dm = Y(rightLeaves,leftLeaves);
if groupCriteria
    srdm=sum(dm,2); scdm=sum(dm,1);
    dm = srdm(:,ones(numLeft,1))+scdm(ones(numRight,1),:);
end
if numLeft == 1
    % if numRight == 1
    % We process the case that the last branches has just two leaves
    % in the beginning. So when we reach hear, numRight must be
    % greater than one.
    
    % We are at a node that has one leaf on the left and a branch
    % on the right; use a vectorized approach to find a right-inner
    % leaf that leads to a maximum score for every possible
    % right-outer leaf:
     [maxR,rIndex] = max((scores(rightLeaves,rightLeaves)),[],1);
     %[~,tR]=max(maxR+dm');% It has similar efficiency as the following
     %loop
     tR = 0;
     best = -inf; 
     for innerR = 1:numRight
         val = dm(innerR) + maxR(innerR);
         if val > best
             tR = innerR;
             best = val;
         end
         
     end
      bestL = leftLeaves;
      right = rIndex(tR);
      bestR =rightLeaves(right);
      tracebackMatrix( bestR,bestL) = rightLeaves(tR);
      tracebackMatrix(bestL,bestR) = bestL;
        
elseif numRight == 1
        % We are at a node that has one leaf on the right and a branch on
        % the left; use a vectorized approach to find a left-inner leaf
        % that leads to a maximum score for every possible left-outer leaf:
       [maxL,LIndex] = max((scores(leftLeaves,leftLeaves)), [],1);
       tL = 0;
       best = -inf;
       for innerL = 1:numLeft
           val = dm(innerL) + maxL(innerL);
           if val > best
                 tL = innerL;
                 best = val;
           end
               
       end
      left = LIndex(tL);
      bestL = leftLeaves(left);
      bestR =rightLeaves;
      tracebackMatrix( bestL,bestR) = leftLeaves(tL);
      tracebackMatrix(bestR,bestL) = bestR;
else %There are more than one nodes on both right and left branches
    [maxL,lIndex] = max((scores(leftLeaves,leftLeaves)),[],1);
    [maxR,rIndex] = max((scores(rightLeaves,rightLeaves)),[],1);
    tL = 0; tR = 0;
    best = -inf;
     for innerL = 1:numLeft
         for innerR= 1:numRight
             val  = maxL(innerL)+dm(innerR,innerL)+maxR(innerR);
             if val > best 
                 best =val;
                 tL = innerL; 
                 tR = innerR;  
             end
               
         end
     end
     left = lIndex(tL);
     right = rIndex(tR);
     bestL = leftLeaves(left); %the best left most leaf
     bestR = rightLeaves(right);%the best right most leaf
     tracebackMatrix(bestR,bestL) = rightLeaves(tR);
     tracebackMatrix(bestL, bestR) = leftLeaves(tL);
end
% Follow the traceback matrix to create a predecessor list 
lq = zeros(numLeaves,1); % left queue
rq = zeros(numLeaves,1); % right queue
dl = zeros(numLeaves,1); % descendants list

i = 1; lq(1) = bestL;  rq(1) = bestR;         % initialize queue
while i>0 % until queue is empty
    left = lq(i); right = rq(i); i = i - 1;   % remove next from queue
    % find inner leaves in the traceback matrix
    lr = tracebackMatrix(left,right); %lr is rightmost in left tree
    rl = tracebackMatrix(right,left); %leftmost in the right tree
    if right~=rl
        i = i + 1; lq(i) = rl; rq(i) = right; % update queue
        dl([lr rl]) = [rl right];             % update descendants list
    end
    if left~=lr
        i = i + 1; lq(i) = left; rq(i) = lr;  % update queue
        dl([left lr]) = [lr rl];              % update descendants list
    end
end

% recover order from the predecessor list
order = zeros(1,numLeaves); 
order(1) = bestL;
for i = 1:numLeaves-1
    order(i+1) = dl(order(i));
end

% flip the tree 
if sum(abs(order-(1:numLeaves)))>sum(abs(order-(numLeaves:-1:1)))
    order = fliplr(order);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nested function  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function kids = getChildrenFromNode(node)
        if node <= numLeaves % nothing to do if node is a leaf
            kids = node;
        else
            % run through the branches of the nodes.
            out = false(1,node);
            out(node) = true;
            for count = node:-1:numLeaves+1
                if out(count)
                    out(Z(count-numLeaves,[1 2])) = true;
                end
            end
            % find the leaves that are children
            kids = find(out(1:numLeaves));
        end
    end % of nested function
end % of main function

