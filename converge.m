function result = converge(func, newres, tol)
    result = newres;
    newres = func(result);
    while (norm(newres - result) > tol)
        result = newres;
        newres = func(result);
    end
end