module ForceAnalysis

using Reexport

using DataFrames
using CategoricalArrays: unique

@reexport using BeForData

export
    # processing
    peak_differences,
    epoch_rejection_ids,
    epoch_rejection,
    aggregate,
    minimum,
    maximum,
    diff,
    # response detection
    OnsetCriterion,
    ForceResponse,
    response_detection,
    peak_force,
    impulse_size,
    duration,
    latency,
    # plotting
    highlight_ranges!,
    plot_av_epoch!,
	plot_good_bad!,
	highlight_ranges!

include("processing.jl")
include("response_detection.jl")

## extensions
_makie_error() = throw(ArgumentError("Have you loaded an appropriate Makie backend?"))
highlight_ranges!(::Any, ::Any, ::Any; kwargs...) = _makie_error()
plot_good_bad!(::Any, ::BeForEpochs; kwargs...) = _makie_error()
plot_av_epoch!(::Any, ::BeForEpochs; kwargs...) = _makie_error()

if !isdefined(Base, :get_extension)
    include("../ext/ForceMakieExt.jl")
end


end
