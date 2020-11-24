# VCF comparison

Ths microservice regroup bash scripts in order to deal with normalisation and comparison of VCF files in quality control and in non regression tests.

## Installation

Use the tool application [docker](https://www.docker.com/) to install the working environment.

```docker
docker-compose up -d


```

## Usage

#TODO
```docker
docker run 
This script is associate with a listener to automate analysis. Thus options will be automatically provide to the script

--run =<RUN>|--cqi=<CQI_SAMPLE>|--genes=<GENES>|--vaf=<VAF,FLOAT,:ALLELE>|--cqigz=<CQI_VCF> [options...]
--run: Full path of run containing STARK results
--cqi: CQI sample among STARK results  
--genes: BED file for the sample  you want to perform analysis. If you don't provide bed it will search one in the sample during analysis. (If the Sample do not contains BED, it will lead to an error)
--vaf: Optionnal option, it will filter VCF according to VAF.  ex --vaf=--min-af,0.10,:major   conf BCFTOOLS view manual)
--cqigz: it represent a list of reference vcf. If any list is present inside the run and the REG option is set to True it will consider the analysis as non regression and a global output will be print in the run. (NB a comparison of all sample according to sample location in the SET of reference will be done)

DEPRECATED:
This microservice is associate with a bash script planned to deal with RTG-tools developp by REALTIME GENOMICS
If the previous script is launched by a listener, then rtg-tools will be launch immediately after. In the case of non regression, you are free to launch this script manually 
Only option --run is usefull. It works as the same as before (full path of the run). It will search for CQI folder in each sample and if one is found it will launch analysis on this one 

#RESULTS
It creates a RES folder inside CQIsample. Inside this folder there are a results file for both SNV and indels stats and a tsv files containing variants miss by STARK pipeline analysis.

NB: It appears that cuting on BED lead to missing some variant. It means that some variants inside REF VCF are located at the edge of BED regions and thereby they are not consider in analysis.
So if you like to consider analysis as the same as STARK V n-1 do not cut on BED to start proces with the same amount of variants

Quality Control

In the sample


```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)

