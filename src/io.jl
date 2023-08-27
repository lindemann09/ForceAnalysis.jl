fl_formats = Dict(:jld2 => ".jld2", :csv => ".csv.zip")

"""
	save(filename::String, fd::ForceData; compress = true)
	save(filename::String, [format::Symbol], fe::ForceEpochs; compress = true)

Saves `ForceData` or `ForceEpochs` as jld2 file.

For `ForceEpochs` the format can be optionally specified:
* `:jld2`: jld2 file
* `:csv`: creates a Zip-file containing seperate csv-data for `forces`, `design`,
		 `baseline` and a json file with `sampling_rate` and `zero_times`
"""
function FileIO.save(filename::String, fd::ForceData; compress = true)
	return jldsave(_check_suffix(filename, :jld2), compress; fd)
end


"""
	load(::Type{ForceData}, filename::String)
	load(::Type{ForceEpochs}, [format::Symbol], filename::String)

Loads `ForceData` or `ForceEpochs` as jdl2 file.

If the format of the `ForceEpochs` is not specified, it will be infered from the
filename suffix.
"""
function FileIO.load(::Type{ForceData}, filename::String)
	return convert(ForceData, load_object(filename))
end;

## Force epochs

function FileIO.save(filename::String, format::Symbol, fe::ForceEpochs; compress = true)
	filename = _check_suffix(filename, format)
	if format == :jld2
		jldsave(filename, compress; fe)
	elseif format == :csv
		save_csv(filename, compress; fe)
	end
end



# jld2 as default save method
function FileIO.save(filename::String, fe::ForceEpochs; compress = true)
	save(filename, :jld2, fe; compress)
end

function FileIO.load(::Type{ForceEpochs}, format::Symbol, filename::String)
	if format == :jld2
		return convert(ForceEpochs, load_object(filename))
	elseif format == :csv
		return load_csv(filename)
	else
		_format_error(format)
	end
end;

# guess format, if not defined
function FileIO.load(::Type{ForceEpochs}, filename::String)
	format = :unknown
	for (k, v) in fl_formats
		if endswith(filename, v)
			format = k
		end
	end
	if format == :unknown
		throw(DomainError(filename, "Please specify file format. Can't guess it from filename."))
	end
	return load(ForceEpochs, format, filename)
end


_csv_error() = throw(ArgumentError(
	"To save CSV, please load the following modules: CSV, JSON and ZipArchives."))
save_csv(::Any, ::ForceEpochs; kwargs...) = _csv_error()
load_csv(::Any; kwargs...) = _csv_error()

## helper
_format_error(format) = throw(DomainError(format,
	"Unkown fileformat. Possible formats are $(keys(fl_formats))"))

function _check_suffix(flname::AbstractString, format::Symbol)
	# ensure correct suffix based on format
	if haskey(fl_formats, format)
		suffix = fl_formats[format]
		return endswith(flname, suffix) ? flname : flname * suffix
	else
		_format_error(format)
	end
end
