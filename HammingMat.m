function HamDistMat = HammingMat(A)
%A is the matrix containing the binary vectors as columns
HamDistMat = zeros(size(A,2),size(A,2));

for i = 1:size(A,2),
    for j = 1:size(A,2),
        HamDistMat(i,j) = sum(abs(A(:,i)-A(:,j)));
    end
end

end