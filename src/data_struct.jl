struct ForceData{T<:FloatOrMissing}
    dat::Vector{T} # force data
    ts::Vector{Int} # timestamps
    sr::Int # sampling_rate
    meta::Dict

    function ForceData(dat::Vector{T}, ts::Vector{Int}, sr::Int, meta::Dict) where {T}
        check_timestamps(ts, length(dat))
        return new{T}(dat, ts, sr, meta)
    end
end;

struct MultiForceData{N,T<:FloatOrMissing}
    dat::Matrix{T} # force data
    ts::Vector{Int} # timestamps
    sr::Int # sampling rate
    ids::NTuple{N,Symbol} # ids/labels
    meta::Dict

    function MultiForceData(dat::Matrix{T}, ts::Vector{Int}, sr::Int,
        ids::Vector, meta::Dict) where {T}
        check_timestamps(ts, size(dat, 1))
        N = size(dat, 2)
        ls = length(ids)
        if ls == 0
            ids = [Symbol(Char(x)) for x in 97:(97 + N - 1)]
        else
            ls == size(dat, 2) || throw(
                ArgumentError(
                    "if ids are defined, they must one label for each sensor."),
            )
            ids = Symbol.(ids)
        end
        return new{N,T}(dat, ts, sr, NTuple{N,Symbol}(ids), meta)
    end
end;

struct ForceProfiles{T<:FloatOrMissing,B<:FloatOrMissing}
    dat::Matrix{T}
    sr::Int # sampling rate
    design::DataFrame
    baseline::Vector{B}
    zero_sample::Int

    function ForceProfiles(force::Matrix{T}, sr::Int, design::DataFrame,
        baseline::Vector{B}, zero_sample::Int) where {T,B}
        lf = size(force, 1)
        lb = length(baseline)
        lf == lb || throw(
            ArgumentError(
                "Number of rows of force ($(lf)) must match the length of baseline ($(lb)).",
            ),
        )
        return new{T,B}(force, sr, design, baseline, zero_sample::Int)
    end
end;

function ForceData(x::MultiForceData, id::Union{Integer,Symbol,String})
    return ForceData(force(x, id), x.ts, x.sr, x.meta)
end

function MultiForceData(
    df::DataFrame;
    sampling_rate::Integer,
    time_stamp_col::Union{Nothing,Symbol,String}=nothing,
    meta::Dict=Dict(),
)
    # dataframe to MultiForceData
    if time_stamp_col isa Symbol
        time_stamp_col = String(time_stamp_col)
    end
    ts = Int[]
    ids = String[]
    for n in names(df)
        if n == time_stamp_col
            ts = df[:, n]
        else
            push!(ids, n)
        end
    end
    return MultiForceData(Matrix(df[:, ids]), ts, sampling_rate, ids, meta)
end

timestamps(x::Union{ForceData,MultiForceData}) = x.ts
sampling_rate(x::Union{ForceProfiles,ForceData,MultiForceData}) = x.sr
n_samples(x::Union{ForceData,MultiForceData}) = size(x.dat, 1)
n_samples(fp::ForceProfiles) = size(fp.dat, 2)
n_profiles(fp::ForceProfiles) = size(fp.dat, 1)

force(x::ForceData) = x.dat
force(x::ForceProfiles) = x.dat
force(x::MultiForceData, id::Int) = x.dat[:, id]
force(x::MultiForceData, id::String) = force(x, Symbol(id))
function force(x::MultiForceData, id::Symbol)
    # returns force
    idx = findfirst(x.ids .== id)
    idx isa Integer || throw(ArgumentError(
        "Can not find force data labelled '$(id)'"))
    return force(x, idx)
end

function duration(n_samples::Int; sampling_rate::Int)
    # calculates n of samples into milliseconds duration
    return n_samples * (1000/sampling_rate)
end


function copy(fp::ForceProfiles)
    return ForceProfiles(copy(fp.dat), fp.sr, copy(fp.design),
        copy(fp.baseline), fp.zero_sample)
end


function Base.show(io::IO, mime::MIME"text/plain", x::ForceData)
    println(io, "ForceData")
    print(io, " $(n_samples(x)) samples, sampling rate: $(sampling_rate(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::MultiForceData)
    println(io, "MultiForceData")
    print(io, "  $(n_samples(x)) samples, sampling rate: $(sampling_rate(x))")
end;


function Base.show(io::IO, mime::MIME"text/plain", x::ForceProfiles)
    println(io, "ForceProfile")
    println(io, "  $(n_profiles(x)) profiles")
    println(io, "  $(n_samples(x)) samples, sampling rate: $(sampling_rate(x)), zero sample: $(x.zero_sample)")
    print(io, "  Design: $(names(x.design))")
end;

# helper
function check_timestamps(timestamps::Vector, required_length::Int)
    lt = length(timestamps)
    return lt == 0 || lt == required_length ||
           throw(
               ArgumentError(
                   "if timestamps are defined, they must have the same length as samples of one sensor.",
               ),
           )
end