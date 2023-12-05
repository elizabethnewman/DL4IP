function[X] = normalizeData(X,normType)

switch normType
    case 'Frobenius'
        X = X / norm(X(:));

    case 'max'
        X = X / max(X(:));

    case 'rescale' % between 0 and 1
        X = rescale(X);
end

end
