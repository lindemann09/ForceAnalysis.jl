
struct ForceExpData
    subject_id::Integer
    file::AbstractString
    force_df::DataFrame
    expdata::DataFrame
    udp::DataFrame
end;

mutable struct ForceProfiles{T<:FloatOrMissing}
    const force::AbstractMatrix{T}
    design ::DataFrame
    const baseline::AbstractVector{T}
    const zero_sample::Integer
end;

copy(fp::ForceProfiles) = ForceProfiles(fp.force, copy(fp.design),
                    fp.baseline, fp.zero_sample)

ncol(fp::ForceProfiles) = size(fp.force, 2)