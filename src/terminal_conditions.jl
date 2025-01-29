function time_termination(c, iteration, f_diff, time)
    if time > c
        return true
    end
    return false
end
function iteration_termination(c, iteration, f_diff, time)
    if iteration > c
        return true
    end
    return false
end
function progress_termination(c, iteration, f_diff, time)
    if f_diff < c
        return true
    end
    return false
end