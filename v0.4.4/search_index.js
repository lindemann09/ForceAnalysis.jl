var documenterSearchIndex = {"docs":
[{"location":"#ForceAnalysis.jl","page":"ForceAnalysis.jl","title":"ForceAnalysis.jl","text":"","category":"section"},{"location":"","page":"ForceAnalysis.jl","title":"ForceAnalysis.jl","text":"","category":"page"},{"location":"","page":"ForceAnalysis.jl","title":"ForceAnalysis.jl","text":"Modules = [ForceAnalysis]","category":"page"},{"location":"#ForceAnalysis.ForceData","page":"ForceAnalysis.jl","title":"ForceAnalysis.ForceData","text":"ForceData{T <: AbstractFloat}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"#ForceAnalysis.ForceEpochs","page":"ForceAnalysis.jl","title":"ForceAnalysis.ForceEpochs","text":"ForceEpochs{T <: AbstractFloat}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"#ForceAnalysis.ForceResponse","page":"ForceAnalysis.jl","title":"ForceAnalysis.ForceResponse","text":"ForceResponse\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"#ForceAnalysis.MultiForceData","page":"ForceAnalysis.jl","title":"ForceAnalysis.MultiForceData","text":"MultiForceData{T <: AbstractFloat}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"#ForceAnalysis.OnsetCriterion","page":"ForceAnalysis.jl","title":"ForceAnalysis.OnsetCriterion","text":"OnsetCriterion{T <: AbstractFloat}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"#Base.copy-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Base.copy","text":"copy(fe::ForceEpochs)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#Base.maximum-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Base.maximum","text":"maximum(fe:ForceEpochs)\n\nMinimum of each epoch.\n\n\n\n\n\n","category":"method"},{"location":"#Base.minimum-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Base.minimum","text":"minimum(fe:ForceEpochs)\n\nMinimum of each epoch.\n\n\n\n\n\n","category":"method"},{"location":"#DataFrames.aggregate-Union{Tuple{ForceEpochs{T}}, Tuple{T}} where T<:AbstractFloat","page":"ForceAnalysis.jl","title":"DataFrames.aggregate","text":"aggregate(fe::ForceEpochs; condition::ColumnIndex = :all,\n\tsubject_id::Union{Nothing, ColumnIndex} = nothing, agg_fnc::Function = mean)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#DataFrames.subset-Tuple{ForceEpochs, Union{Tuple{Vararg{Integer}}, AbstractVector{<:Integer}}}","page":"ForceAnalysis.jl","title":"DataFrames.subset","text":"subset(fe::ForceEpochs, rows::Base.AbstractVecOrTuple{Integer})\nsubset(fe::ForceEpochs, args...)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#FileIO.load-Tuple{Type{ForceData}, String}","page":"ForceAnalysis.jl","title":"FileIO.load","text":"load(::Type{ForceData}, filename::String)\nload(::Type{ForceEpochs}, [format::Symbol], filename::String)\n\nLoads ForceData or ForceEpochs as jdl2 file.\n\nIf the format of the ForceEpochs is not specified, it will be infered from the filename suffix.\n\n\n\n\n\n","category":"method"},{"location":"#FileIO.save-Tuple{String, ForceData}","page":"ForceAnalysis.jl","title":"FileIO.save","text":"save(filename::String, fd::ForceData; compress = true)\nsave(filename::String, [format::Symbol], fe::ForceEpochs; compress = true)\n\nSaves ForceData or ForceEpochs as jld2 file.\n\nFor ForceEpochs the format can be optionally specified:\n\n:jld2: jld2 file\n:csv: creates a Zip-file containing seperate csv-data for forces, design,\n\n\t `baseline` and a json file with `sampling_rate` and `zero_times`\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.adjust_baseline!-Tuple{ForceEpochs, UnitRange{<:Integer}}","page":"ForceAnalysis.jl","title":"ForceAnalysis.adjust_baseline!","text":"adjust_baseline!(fe::ForceEpochs, baseline_window::UnitRange{<:Integer})\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.concatenate-Tuple{ForceEpochs, ForceEpochs}","page":"ForceAnalysis.jl","title":"ForceAnalysis.concatenate","text":"concatenate(a::ForceEpochs, b::ForceEpochs)\n\nconcatenate two or multiple ForceEpochs\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.epoch_rejection-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"ForceAnalysis.epoch_rejection","text":"epoch_rejection(fe::ForceEpochs; force_range::UnitRange, max_difference::Integer,\n\tmax_diff_windows_size::Integer)\n\nReturn ForceEpochs fitting the criteria for 'good' epochs\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.epoch_rejection_ids-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"ForceAnalysis.epoch_rejection_ids","text":"epoch_rejection_ids(fe::ForceEpochs; force_range::UnitRange, max_difference::Integer,\n\tmax_diff_windows_size::Integer)\n\nReturn BitVector indicating 'good' epochs\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.epochs-Union{Tuple{ForceData{T}}, Tuple{T}} where T<:AbstractFloat","page":"ForceAnalysis.jl","title":"ForceAnalysis.epochs","text":"epochs(fd::ForceData{<: AbstractFloat};\n\tzero_times::AbstractVector{<:Integer},\n\tn_samples::Integer,\tn_samples_before::Integer)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.force-Tuple{ForceData}","page":"ForceAnalysis.jl","title":"ForceAnalysis.force","text":"force(fd::ForceData)\nforce(fd::MultiForceData, id::Union{Int, String, Symbol})\nforce(fe::ForceEpochs)\n\nReturns the force data as Matrix or Vector.\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.impulse_size-Tuple{AbstractVector{<:AbstractFloat}, ForceResponse}","page":"ForceAnalysis.jl","title":"ForceAnalysis.impulse_size","text":"impulse_size(force_vector::AbstractVector{<:AbstractFloat}, rb::ForceResponse)\nimpulse_size(fp::ForceEpochs, rb::AbstractVector{ForceResponse})\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.latency-Tuple{ForceResponse}","page":"ForceAnalysis.jl","title":"ForceAnalysis.latency","text":"latency(force_vector::AbstractVector{<:AbstractFloat}, rb::ForceResponse)\nlatency(fp::ForceEpochs, rb::AbstractVector{ForceResponse})\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.lowpass_filter!-Tuple{ForceData}","page":"ForceAnalysis.jl","title":"ForceAnalysis.lowpass_filter!","text":"lowpass_filter!(fd::ForceData; cutoff_freq::Real,butterworth_order::Integer = 4)\nlowpass_filter!(fd::MultiForceData; cutoff_freq::Real,butterworth_order::Integer = 4)\nlowpass_filter!(fd::ForceEpochs; cutoff_freq::Real,butterworth_order::Integer = 4)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.response_detection-Tuple{AbstractVector{<:AbstractFloat}, OnsetCriterion}","page":"ForceAnalysis.jl","title":"ForceAnalysis.response_detection","text":"response_detection(force_vector::AbstractVector{<:AbstractFloat}, criterion::OnsetCriterion; zero_sample::Integer)::ForceResponse\nresponse_detection(fe::ForceEpochs, criterion::OnsetCriterion)::Vector{ForceResponse}\n\nReturns a ForceResponse or Vector{ForceResponse} TODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.scale_force!-Tuple{ForceData, Real}","page":"ForceAnalysis.jl","title":"ForceAnalysis.scale_force!","text":"scale_force!(fd::ForceData, factor::Real)\nscale_force!(fe::ForceEpochs, factor::Real)\n\nTODO\n\n\n\n\n\n","category":"method"},{"location":"#ForceAnalysis.sderr-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"ForceAnalysis.sderr","text":"sderr(fe::ForceEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)\n\nComputes the standard error of the epochs, i.e., column-wise standard deviations for each sample. If rows is defined, only the selected rows are considered.\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.mean-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Statistics.mean","text":"mean(fe::ForceEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)\n\nComputes the mean of the epochs, i.e., column-wise means for each sample. If rows is defined, only the selected rows are considered.\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.median-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Statistics.median","text":"median(fe::ForceEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)\n\nComputes the median of the epochs, i.e., column-wise medians for each sample. If rows is defined, only the selected rows are considered.\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.std-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Statistics.std","text":"std(fe::ForceEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)\n\nComputes the standard deviation of the epochs, i.e., column-wise standard deviations for each sample. If rows is defined, only the selected rows are considered.\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.var-Tuple{ForceEpochs}","page":"ForceAnalysis.jl","title":"Statistics.var","text":"var(fe::ForceEpochs; rows::Union{Nothing, BitVector, Vector{<:Integer}}=nothing)\n\nComputes the variance of the epochs, i.e., column-wise variances for each sample. If rows is defined, only the selected rows are considered.\n\n\n\n\n\n","category":"method"}]
}