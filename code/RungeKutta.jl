function rk4(f, init_condition, init_time, stop_condition; params=nothing, step_size = 0.01)
    t = init_time
    y = init_condition
    ts = [t]
    ys = [y]
    # log_point = 0.2
    while !stop_condition(t, y)
        k1 = f(t, y, params)
        # vector operations on y
        k2 = f(t + step_size/2, y .+ step_size .* k1 ./2, params) 
        k3 = f(t + step_size/2, y .+ step_size .* k2 ./2, params)
        k4 = f(t + step_size, y .+ step_size .* k3, params)
        @. y += step_size/6 * (k1 + k2 * 2 + k3 * 2 + k4) # vector ops on y
        push!(ys, y) # append new y to the history of y
        t += step_size
        push!(ts, t) # append new t to the history of t. Good to have in case step size is adaptive.
        # if y[1] > log_point
        #     println((t,y))
        #     log_point += 0.1
        # end
    end
    return (ts, ys)
end


function rk4_update(f, init_condition, init_time, stop_condition; update_size=10, step_size = 0.01)
    t = init_time
    y = init_condition
    ch = Channel(update_size) do ch
        while !stop_condition(t, y)
            k1 = f(t, y)
            # vector operations on y
            k2 = f(t + step_size/2, @. y + step_size * k1/2) 
            k3 = f(t + step_size/2, @. y + step_size * k2/2)
            k4 = f(t + step_size, @. y + step_size * k3)
            @. y += step_size/6 * (k1 + k2 * 2 + k3 * 2 + k4) # vector ops on y
            t += step_size
            put!(ch, (t,y))
        end
    end
    return ch
end