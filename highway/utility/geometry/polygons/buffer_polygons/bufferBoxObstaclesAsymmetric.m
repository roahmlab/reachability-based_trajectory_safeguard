function Bout = bufferBoxObstaclesAsymmetric(B,L,W)
% Given a 2-by-n matrix of 2-D box obstacles, extend them by a distance L
% on their longer dimension and W on their shorter dimension

    N = 13 ;
    
    % check beginning and end of B for nans
    if ~isnan(B(1,end))
        B = [B, nan(2,1)] ;
    end
    
    if isnan(B(1,1))
        B = B(:,2:end) ;
    end

    % not gonna bother making circles on the corners, just stretch the
    % boxes out boiiii
    Bout = nan(size(B)) ;
    
    % for each box...
    for idx = 1:6:size(B,2)-1
        % get the box
        Bidx = B(:,idx:idx+5) ;
        
        % find the box's side widths
        b1 = Bidx(:,1) ; b2 = Bidx(:,2) ; b3 = Bidx(:,3) ;
        Ltemp = sqrt(sum((b2 - b1).^2)) ;
        Wtemp = sqrt(sum((b3 - b2).^2)) ;
        
        % make sure we've got the dimensions right
        if Ltemp > Wtemp
            Lidx = Ltemp ;
            Widx = Wtemp ;
        else
            Lidx = Wtemp ;
            Widx = Ltemp ;
        end
        
        % find the box's centroid and rotation angle
        c = mean(Bidx(:,1:4),2) ;
        h = atan2(b2(2) - b1(2),b2(1) - b1(1))  ;
        
        % make a new box with the same centroid and angle but the desired
        % length and width
        Bnew = rotmat(h)*makeBox([Lidx + 2*L, Widx + 2*W]) + repmat(c,1,5) ;
        Bout(:,idx:idx+5) = [Bnew, nan(2,1)] ;
    end
    
    % remove final nan column
    Bout = Bout(:,1:end-1) ;
end