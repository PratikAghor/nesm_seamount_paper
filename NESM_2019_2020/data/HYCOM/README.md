* generate directory for each field first, for example, ```mkdir salinity```
* ```sbatch gethycom_salinity.sbatch``` will download all the data for each day in the salinity directory
* Repeat for all the other fields
* use ```../merged/merge_hycom_data.sbatch``` to combine individual files into one .nc file
