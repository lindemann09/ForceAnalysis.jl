
struct ForceExpData
    subject_id::Integer
    file::AbstractString
    force_df::DataFrame
    expdata::DataFrame
    udp::DataFrame
end;

mutable struct ForceProfiles{T<:FloatOrMissing, B<:FloatOrMissing}
    const force::AbstractMatrix{T}
    design ::DataFrame
    const baseline::AbstractVector{B}
    const zero_sample::Integer
end;

copy(fp::ForceProfiles) = ForceProfiles(copy(fp.force), copy(fp.design),
                    copy(fp.baseline), fp.zero_sample)

ncol(fp::ForceProfiles) = size(fp.force, 2)