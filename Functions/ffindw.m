function  [w_index, singular] = ffindw(A, ratio)
w_index = 0;
singular = svd(A);
for i = 1 : rank(A)
    if sum(singular(1 : i)) > ratio * sum(singular)
        w_index = i;
        break;
    end
end
end