module ForceAnalysis

using UnPack
using DSP # signal processing, filtering
using DataFrames
using JLD2, CodecZlib
using FileIO
using CategoricalArrays: unique

using Statistics

export ForceData, # force data
    MultiForceData,
    ForceEpochs,
    force,
    # io
    save,
    load,
    # preprocessing
    scale_force!,
    lowpass_filter!,
    epochs,
    adjust_baseline!,
    # processing
    peak_differences,
    epoch_rejection,
    aggregate,
    subset,
    # response detection
    OnsetCriterion,
    ForceResponse,
    response_detection,
    peak_force,
    impulse_size,
    duration,
    latency,
    # force statistics
    minimum,
    maximum,
    mean,
    std,
    var,
    median,
    diff,
    # plotting
    highlight_ranges!,
    plot_av_epoch!,
	plot_good_bad!,
	highlight_ranges!

include("data_structs.jl")
include("io.jl")
include("stats.jl")
include("preprocessing.jl")
include("processing.jl")
include("response_detection.jl")

## extensions
_makie_error() = throw(ArgumentError("Have you loaded an appropriate Makie backend?"))
highlight_ranges!(::Any, ::Any, ::Any; kwargs...) = _makie_error()
plot_good_bad!(::Any, ::ForceEpochs; kwargs...) = _makie_error()
plot_av_epoch!(::Any, ::ForceEpochs; kwargs...) = _makie_error()

if !isdefined(Base, :get_extension)
    include("../ext/ForceMakieExt.jl")
end

if !isdefined(Base, :get_extension)
    include("../ext/CSVFilesExt.jl")
end


end
