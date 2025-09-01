function Base.diff(fe::BeForEpochs; dims::Integer)
    mtx = fe.dat
    if dims == 1
        z = zeros(Float64, 1, size(mtx, 2))
        dat = vcat(z, diff(mtx; dims=1))
    elseif dims == 2
        z = zeros(Float64, size(mtx, 1), 1)
        dat = hcat(z, diff(mtx; dims=2))
    else
        throw(ArgumentError("dims has to be 1 or 2 and not $dims"))
    end
    return BeForEpochs(dat, fe.sampling_rate, fe.design, fe.baseline,
            fe.zero_sample,  meta=copy(fe.meta))
end;


"""
    minimum(fe:BeForEpochs)

Minimum of each epoch.
"""
Base.minimum(fe::BeForEpochs) = return vec(minimum(fe.dat, dims=2))

"""
    maximum(fe:BeForEpochs)

Minimum of each epoch.
"""
Base.maximum(fe::BeForEpochs) = return vec(maximum(fe.dat, dims=2))
