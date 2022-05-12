function munge(sol)
    M = zeros(length(sol),length(sol[1]))
    for i = 1:length(sol)
        M[i,:] = sol[i]
    end
    return M
end
