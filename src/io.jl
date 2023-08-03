fl_formats = Dict(:jld2 => ".jld2", :csv => ".csv.zip")

function save(filename::String, forces::ForceData)
	return jldsave(_check_suffix(filename, :jld2), true; forces)
end

function load(::Type{ForceData}, filename::String)
	return convert(ForceData, load_object(filename))
end;

## Force epochs

function save(filename::String, format::Symbol, fe::ForceEpochs; compress = true)
	filename = _check_suffix(filename, format)
	if format == :jld2
		jldsave(filename, compress; fe)

	elseif format == :csv
		dataname = _dataname(filename)
		ZipWriter(filename) do zipfl
			zip_newfile(zipfl, dataname * ".forces.csv"; compress)
			CSV.write(zipfl, CSV.Tables.table(fe.dat); writeheader = false)

			zip_newfile(zipfl, dataname * ".design.csv"; compress)
			CSV.write(zipfl, fe.design)

			zip_newfile(zipfl, dataname * ".baseline.csv"; compress)
			CSV.write(zipfl, CSV.Tables.table(fe.baseline); writeheader = false)

			zip_newfile(zipfl, dataname * ".meta.json"; compress)
			JSON.print(zipfl, Dict(:sampling_rate => fe.sampling_rate,
					:zero_sample => fe.zero_sample), 1)
		end
	end
	return nothing
end

# jld2 as default save method
function save(filename::String, fe::ForceEpochs; compress = true)
	save(filename, :jld2, fe; compress)
end

function load(::Type{ForceEpochs}, format::Symbol, filename::String)
	if format == :jld2
		return convert(ForceEpochs, load_object(filename))
	elseif format == :csv
		dataname = _dataname(filename)

		zipfl = zip_open_filereader(filename)
        txt = zip_readentry(zipfl, dataname * ".meta.json", String)
        meta = JSON.parse(txt)
        fl = zip_openentry(zipfl, dataname * ".forces.csv")
        dat = CSV.read(fl, CSV.Tables.matrix; header = false)
        fl = zip_openentry(zipfl, dataname * ".baseline.csv")
        baseline = CSV.read(fl, CSV.Tables.matrix; header = false)
        fl = zip_openentry(zipfl, dataname * ".design.csv")
        design = CSV.read(fl, DataFrame)
        close(zipfl)

        return ForceEpochs(dat, meta["sampling_rate"], design,
			vec(baseline), meta["zero_sample"])
	else
		_format_error(format)
	end
end;

# guess format, if not defined
function load(::Type{ForceEpochs}, filename::String)
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


## helper
_format_error(format) = throw(DomainError(format,
	"Unkown fileformat. Possible formats are $(keys(fl_formats))"))

_dataname(flname) = first(split(last(splitdir(flname)), "."))

function _check_suffix(flname::AbstractString, format::Symbol)
	# ensure correct suffix based on format
	if haskey(fl_formats, format)
		suffix = fl_formats[format]
		return endswith(flname, suffix) ? flname : flname * suffix
	else
		_format_error(format)
	end
end
