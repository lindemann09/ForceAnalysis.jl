module ForceAnalysis

using Reexport

using DataFrames
using CategoricalArrays: unique

@reexport using BeForData

export
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
