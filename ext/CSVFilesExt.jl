module CSVFilesExt

using CSV, JSON
using ZipArchives
using DataFrames

using ForceAnalysis

function ForceAnalysis.save_csv(filename::String, fe::ForceEpochs; compress = true)
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
	return nothing
end

function ForceAnalysis.load_csv(filename::String)
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
end

_dataname(flname) = first(split(last(splitdir(flname)), "."))

end