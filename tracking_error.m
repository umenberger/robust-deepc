function te = tracking_error(y,r,n)

    e = y - r(:,1:size(y,2));
    
    te = norm(e(:),n); 
    
end
    