function [Q, R] = mGS(A)
[m,n] = size(A);

if m < n
    error('matrix dimension ineligible (m<n)')
end
Q = zeros(m,n); R = zeros(n,n);

R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)./R(1,1);

for j = 2:n    
    R(1:(j-1), j) = conj(Q(:, 1:(j-1))') * A(:,j);
    qj = A(:,j) - Q(:, 1:(j-1)) * R(1:(j-1), j);
    R(j,j) = norm(qj);
    Q(:,j) = qj ./ R(j,j);
    
    %{
    if mod(j,100)==0
        disp(j)
    end
    %}
end

return