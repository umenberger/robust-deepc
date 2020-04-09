function H = vec2block_hankel(x,d)

[nx,L] = size(x);

nc = L - d + 1;

H = zeros(nx*d,nc);

for i = 1:nc
    
    tmp = x(:,(i-1)+1:i+d-1);
    
    H(:,i) = tmp(:);
    
end

