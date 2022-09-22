from operator import index
import pandas as pd
from tqdm import tqdm
import os
import subprocess
from multiprocessing import Pool
from os.path import join as osj

config = {
    "tools": {
        "java": "java",
        "gatk": "/STARK/tools/gatk/3.8-1-0/bin/GenomeAnalysisTK.jar",
        "VaRank": "/home1/TOOLS/tools/varank/VaRank_1.4.3",
        "Alamut-batch": "/home1/TOOLS/tools/alamut_batch/alamut-batch-standalone-1.11",
        "howard": "/STARK/tools/howard/current/bin/HOWARD",
        "bcftools": "/STARK/tools/bcftools/current/bin/bcftools",
        "bgzip": "/STARK/tools/htslib/current/bin/bgzip",
        "tabix": "/STARK/tools/htslib/current/bin/tabix",
        "threads": 2,
    },
    "family": {
        "FATHER": {
            "sex": "M",
            "affinity": ["ASG2104494", "FATHER"],
            "bam": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/ASG2104494/STARK/ASG2104494.bwamem.bam",
            "vcf": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/ASG2104494/STARK/ASG2104494.reports/ASG2104494.final.vcf.gz",
        },
        "MOTHER": {
            "sex": "F",
            "affinity": ["ASG2104493", "MOTHER"],
            "bam": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/ASG2104493/STARK/ASG2104493.bwamem.bam",
            "vcf": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/ASG2104493/STARK/ASG2104493.reports/ASG2104493.final.vcf.gz",
        },
    },
    "foetus": {
        "name": "FAV2104492",
        "bam": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/FAV2104492/STARK/FAV2104492.bwamem.bam",
        "vcf": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/FAV2104492/STARK/FAV2104492.reports/FAV2104492.final.vcf.gz",
    },
    "env": {
        "output": "/STARK/output/res",
        "genome": "/STARK/databases/genomes/current/hg19.fa",
        "bed": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/ASG2104494/STARK/ASG2104494.bed",
        "run": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/",
        "repository": "/STARK/data/users/jb/210728_NB551027_0901_AHLG3NBGXG/",
        "depository": "",
    },
}


def test(result):
    print("#[INFO] DENOVOFILTER")
    # print(result)
    print(*result.columns)
    result["denovo_filter"] = None
    # result["Base_counts_FATHER"] = result["Base_counts_FATHER"].astype(str)
    for i, var in result.iterrows():
        ref = var["REF"]
        mut = var["ALT"]
        freq = [
            item.split(":") for item in var["Base_counts_FATHER"].split("|")[0].split()
        ]
        # print(result["Base_counts_FATHER"].dtypes)
        # print(var["Base_counts_FATHER"].split("|"))
        # for item in var["Base_counts_FATHER"].split("|")[0].split():
        #    print(item)
        #    print(type(item))
        # exit()
        # print(freq)
        freq.extend(
            [
                item.split(":")
                for item in var["Base_counts_MOTHER"].split("|")[0].split()
            ]
        )
        # [dict(item.split(':')[0] = item.split(':')[1]) for item in freq_f]
        val = []
        # print(freq)
        for values in freq:
            # print(values)
            # time.sleep(10)
            if values[0] == mut:
                # print(values[0], mut)
                val.append(int(values[1]))
        # filter only for SNV
        if (
            val
            and all(value < 3 for value in val)
            and var["Total_Depth_FATHER"] >= 300
            and var["Total_Depth_MOTHER"] >= 300
            and len(mut) == 1
            and len(ref) == 1
        ):
            result.loc[i, "denovo_filter"] = "Denovo"
            # print("#[INFO] counts", val)
            # print(result.loc[i, 'denovo_filter'])
            # time.sleep(5)
        elif (
            val
            and all(value < 2 for value in val)
            and var["Total_Depth_FATHER"] >= 200
            and var["Total_Depth_MOTHER"] >= 200
            and len(mut) == 1
            and len(ref) == 1
        ):
            result.loc[i, "denovo_filter"] = "LikelyDenovo"
            # print("#[INFO] counts", val)
            # print(result.loc[i, 'denovo_filter'])
            # print(var)
        elif len(mut) > 1 and len(ref) == 1:
            result.loc[i, "denovo_filter"] = "InDel"
        elif len(ref) > 1 and len(mut) == 1:
            result.loc[i, "denovo_filter"] = "InDel"
        else:
            result.loc[i, "denovo_filter"] = "Inherit"
        # print("#LINE ", [var["#CHROM"], var["POS"], var["REF"], var["ALT"]])
    print(result["denovo_filter"].value_counts())
    return result


def uniqueChr(df):
    unique = df["#CHROM"].unique()
    return unique


def addFrequence(df, sample, s_list):
    # df.sort_values(by=['#CHROM', 'POS'], ascending=[False, True], key=lambda x: np.argsort(index_natsorted(df['#CHROM']
    locus = "Locus_" + sample
    data = []
    # print(df.columns)
    # print(df.head())
    for i, values in df.iterrows():
        # print(values)
        last = {}
        last["#CHROM"] = values[locus].strip().split(":")[0]
        last["POSITION"] = values[locus].strip().split(":")[1]
        # for each indivi stats
        for indiv in s_list:
            tmp = values["Base_counts_" + indiv].strip().split()
            tmp_data = {}
            tmp2_data = {}
            # details of each bases count
            for bases in tmp:
                base = bases.split(":")
                tmp_data[base[0]] = base[1]
                # print(tmp_data)
            for count in tmp_data.items():
                # print(count)
                try:
                    tmp2_data[count[0]] = str(
                        round(float(count[1]) / values["Total_Depth_" + indiv], 2)
                    )
                except:
                    # print("#[INFO] Depth is 0 for ", values["Locus_" + indiv])
                    tmp2_data[count[0]] = "0"

            count = "|".join(":".join((key, val)) for (key, val) in tmp_data.items())
            freq = "|".join(":".join((key, val)) for (key, val) in tmp2_data.items())

            last["Depth_" + indiv] = values["Depth_" + indiv]
            last["Base_counts_" + indiv] = count
            last["Base_frequence_" + indiv] = freq

        data.append(last)
    final = pd.DataFrame(data)
    return final


def createFinaldf(df, header, oup):
    raw = df[df.columns[0:10]]
    with open(oup, "w+") as f:
        for lines in header:
            f.write(lines.strip() + "\n")
    raw.to_csv(oup, mode="a", index=False, header=False, sep="\t")
    return raw


def parse_sample_field(dfVar):
    # UPDATE 21/06/2022
    # handle multisample and keep information of sample field like <sample>_<field> for each annotations

    #############
    ### Parsing Sample Field in VCF
    #############

    dico = []
    # dfTest = pd.Series(dfVar.TEM195660.values,index=dfVar.FORMAT).to_dict()
    bad_annotation = []
    sample_list = []

    # Parsing FORMAT field in VCF
    # print("[#INFO] Parsing FORMAT field")
    isample = list(dfVar.columns).index("FORMAT") + 1
    # index: line where the caller identify an event somethings
    for col in dfVar.columns[isample:]:
        # print("#[INFO] " + col + "\n")
        tmp_ = []
        sample_list.append(col)
        for i, row in tqdm(
            dfVar.iterrows(),
            total=dfVar.shape[0],
            leave=True,
            desc="#[INFO] SAMPLE " + col,
        ):
            if len(row["FORMAT"].split(":")) != len(row[col].split(":")):
                bad_annotation.append(pd.Series(row[:], index=dfVar.columns))
                continue
            else:
                # If more than on sample
                if len(dfVar.columns[isample:]) > 1:
                    index = [col + "_" + x for x in row["FORMAT"].split(":")]
                else:
                    index = row["FORMAT"].split(":")
                toadd = pd.Series(
                    row[col].split(":"),
                    index=index,
                ).to_dict()
                # toadd.update({"caller":col})
                tmp_.append(toadd)
        dico.append(pd.DataFrame(tmp_))

    dfSample = pd.concat(dico, axis=1)
    df_bad_anno = pd.DataFrame(bad_annotation)
    try:
        df_final = dfVar.join(dfSample, how="inner")
    except ValueError:
        df_final = dfVar.join(dfSample, how="inner", lsuffix="_INFO", rsuffix="_SAMPLE")
        print(
            "WARNING some columns are present in both INFO and SAMPLE field, add _INFO and _SAMPLE suffix"
        )
    df_final.drop(columns=sample_list, inplace=True)
    return df_final, df_bad_anno, sample_list


def parseSpeed(covfiles, df, col, sample):
    """
    input: folder containing cov depth file by chr, vcf of the sample and list of unique chr
    """
    # hide warning

    pd.set_option("mode.chained_assignment", None)

    if os.path.basename(covfiles) == "COVFILEchr3":
        print("#[INFO] chr " + covfiles)
    # print("#[INFO] Chrcov "+covfiles)
    data = []
    # read cov file results from gatk depth coverage to dataframe
    tmp = pd.read_csv(covfiles, sep="\t", names=col)
    # chr number
    chr = os.path.basename(covfiles)[7:]
    locus_name = "Locus_" + sample
    cov = tmp.loc[tmp[locus_name].str.startswith(chr + ":")]  # .to_dict    ('records')
    # print("INFO COV", cov.head())
    # select only one chromosome (cut by starmap func to multiprocess)
    df_tmp = df.loc[df["#CHROM"] == chr]
    # to match depth coverage informations
    df_tmp.loc[:, locus_name] = (
        df.loc[:, "#CHROM"].astype(str) + ":" + df.loc[:, "POS"].astype(str)
    )
    # Locus column as index for both variant and depthcoverage file
    df_tmp = df_tmp.set_index(locus_name)
    cov = cov.set_index(locus_name)
    # print("INFO DF_TMP "+covfiles, df_tmp)
    # print("INFO length dftmp ", len(df_tmp.index))
    # Merge variants df and depthcov df in index Locus for each variants    position weget the values of cov
    print("\n")
    print("#[INFO] df tmp", df_tmp)
    print("#[INFO] cov", cov)
    print("\n")
    final = df_tmp.merge(cov, left_index=True, right_index=True, how="left")
    # reset index and reinsert Locus columns ex chr1:5667097 at good column     position
    print("final before index", final)
    final.reset_index(level=0, inplace=True)
    locus = final[locus_name]
    final.drop(columns=locus_name, inplace=True)
    final.insert(10, locus_name, locus)
    final.drop_duplicates(inplace=True)
    final = final.iloc[:, 1:]
    print(*final.columns)
    print(final)
    print("#[INFO] Length variants " + os.path.basename(covfiles), len(final.index))
    return final


def vcfTodataframe(file, rheader=False):
    """
    Take in input vcf file, or tsv and return a dataframe
    I"m gonna build my own vcf parser et puis c'est tout
    return 3 Dataframe, full, only sample, only info
    """
    name, extension = os.path.splitext(file)

    header = []
    variants_tmp = []
    variants = []

    print("#[INFO] " + file)
    if extension == ".vcf":
        # print('[#INFO] VCF: '+file)
        with open(file) as f:
            for lines in f:
                if lines.startswith("#"):
                    header.append(lines.strip())
                else:
                    variants_tmp.append(lines)

    col = header[-1].strip().split("\t")
    for v in variants_tmp:
        variants.append(v.strip().split("\t"))

    rows = []
    # Creating Dataframe from the whole VCF
    print("#[INFO] Whole VCF to Dataframe")
    for i, var in enumerate(variants):
        rows.append(pd.Series(var, index=col))
    # .T necessary to transpose index and columns
    dfVar = pd.concat(rows, axis=1, ignore_index=True).T

    if rheader:
        return dfVar, header
    else:
        return dfVar


def barcode_DPNI(df, foetus):
    dico = {}
    denovo = []
    dico["MOTHER"] = []
    dico["FATHER"] = []
    dico["FF5"] = []
    dico["DPNI_Barcode_coment"] = []
    for i, rows in tqdm(
        df.iterrows(), total=df.shape[0], leave=True, desc="#[INFO] Barcode DPNI ..."
    ):
        # print("INFO index", i)
        for fam in ["MOTHER", "FATHER", "FF5"]:
            tmp = []
            # control var
            if all(
                elem in df.columns
                for elem in [fam + "_GT", fam + "_AD", fam + "_DP", fam + "_FT"]
            ):
                # Variant call heterozygous
                if (
                    rows[fam + "_GT"] == "0/1"
                    or rows[fam + "_GT"] == "1/0"
                    or rows[fam + "_GT"] == "1|0"
                    or rows[fam + "_GT"] == "0|1"
                ):
                    # ensure truethness of var, and no filter in FT (low GQ etc) same for homo
                    if int(rows[fam + "_DP"]) > 30 and (
                        rows[fam + "_FT"] == "."
                        or rows[fam + "_FT"] == "PASS"
                        or rows[fam + "_FT"] == ""
                    ):
                        dico[fam].append("1")
                    else:
                        dico[fam].append("0")
                        tmp.append(
                            "|".join(
                                [
                                    str(rows[fam + "_FT"]),
                                    "Depth:" + str(rows[fam + "_DP"]),
                                ]
                            )
                        )
                # homozygous
                elif rows[fam + "_GT"] == "1/1" or rows[fam + "_GT"] == "1|1":
                    if int(rows[fam + "_DP"]) > 30 and (
                        rows[fam + "_FT"] == "."
                        or rows[fam + "_FT"] == "PASS"
                        or rows[fam + "_FT"] == ""
                    ):
                        dico[fam].append("2")
                    else:
                        dico[fam].append("0")
                        tmp.append(
                            ",".join(
                                [
                                    str(rows[fam + "_FT"]),
                                    "Depth" + fam + ":" + str(rows[fam + "_DP"]),
                                ]
                            )
                        )
                # Ref allele
                else:
                    dico[fam].append("0")
            else:
                print(
                    "ERROR miss one metrics columns in EXIT"
                    + " ".join([fam + "_GT", fam + "_AD", fam + "_DP", fam + "_FT"])
                )
                exit()
        # print(dico)
        # Correct Foetus var
        if (
            dico["MOTHER"][i] == "2"
            and dico["FATHER"][i] == "2"
            and dico["FF5"][i] != "2"
        ):
            tmp.append("GT-issues")
        elif (
            dico["MOTHER"][i] == "0"
            and dico["FATHER"][i] == "0"
            and dico["FF5"][i] != "0"
        ):
            tmp.append("perhaps-denovo")
        else:
            tmp.append("Inherit")
        dico["DPNI_Barcode_coment"].append(tmp)
        # DENOVO filter check as each line
        if (
            str(rows["FF5_AD"]) != "."
            and str(rows["FF5_AD"]) != ""
            and str(rows["FF5_AD"]) != "nan"
            and len(str(rows["FF5_AD"]).split(",")) == 2
        ):
            if int(rows["FF5_AD"].split(",")[1]) < 3 and int(rows["FF5_DP"]) >= 300:
                denovo.append("Denovo")
            elif int(rows["FF5_AD"].split(",")[1]) < 3 and int(rows["FF5_DP"]) >= 200:
                denovo.append("LikelyDenovo")
            else:
                denovo.append("Inherit")
        else:
            denovo.append("Inherit")
    values = [
        dico["MOTHER"],
        dico["FF5"],
        dico["FATHER"],
    ]

    df["DPNI_Barcode"] = ["'" + "".join(item) + "'" for item in list(zip(*values))]
    df["DPNI_Barcode_coment"] = [",".join(item) for item in dico["DPNI_Barcode_coment"]]
    df["DPNI_Denovo"] = denovo
    return df


def add_info(vcf_in, vcf_out, cov):
    sample_name = os.path.basename(vcf_in).split(".")[0]
    if sample_name == config["foetus"]["name"]:
        sample_name = "foetus"
    # transform raw vcf into dataframe and insert header in list
    df, head = vcfTodataframe(vcf_in, True)

    # unique chr cov file from gatk
    unique = uniqueChr(df)
    list_var = []

    # columns of cov file
    with open(cov) as f:
        col = f.readline().strip().split("\t")
    print("#[INFO] Columns depthcov ", col)
    # col = ['Locus', 'Total_Depth', 'Average_Depth', 'Depth', 'Base_counts']
    # parse variants by chr to speed up analysis
    # prepare args for multiprocessing
    print("#[INFO] Parse variants by chr")
    list_cov = [
        (osj(vcf_out + ".dict", covfiles), df, col, sample_name)
        for covfiles in os.listdir(vcf_out + ".dict")
        if not covfiles.startswith("COVFILELocus")
    ]
    # print(list_cov)
    # print("#[INFO] LISTCV", list_cov)
    print("#[INFO] chr nbr", len(list_cov))
    with Pool(6) as pool:
        # TEST async
        result = tqdm(pool.starmap(parseSpeed, list_cov), total=len(list_cov))
        if result.ready():
            print("#[INFO] Starmap done")
            pool.close()
            pool.join()
        for locus in result:
            print(locus.head())
            print(len(locus.index))
            list_var.append(locus)
    # print("#[INFO] LISTVAR ", list_var)
    # result = pd.concat(list_var, axis=1, ignore_index=True)
    print("#[INFO] Process annotations")
    # frequence annotations
    metrics = pd.concat(list_var)
    metrics = metrics.dropna(subset=["REF", "ALT"], how="all")
    # calculate frequence regarding depth and base identification
    s_list = ["foetus"]
    s_list.extend(config["family"])
    print("#[INFO] Family people ", s_list)

    frequence = addFrequence(metrics, sample_name, s_list)
    # remove duplicate column used to create frequence
    metrics.drop(
        columns=[
            "Depth_FATHER",
            "Depth_MOTHER",
            "Depth_foetus",
            "Base_counts_FATHER",
            "Base_counts_MOTHER",
            "Base_counts_foetus",
        ],
        inplace=True,
    )
    dfinal = pd.concat(
        [metrics.reset_index(drop=True), frequence.reset_index(drop=True)], axis=1
    )
    info = frequence.iloc[:, 2:]
    for items in list(info.columns):
        info[items] = items + "=" + info[items].astype(str)
    info["INFOS"] = info.agg(";".join, axis=1)
    dfinal["INFO"] = dfinal["INFO"] + ";" + info["INFOS"]
    # denovofilter
    # dfinal.to_csv("/app/res/test.csv", header=True, sep="\t", index=False)
    # The global denovo filter is more precise, can remove 2 lines above
    # df_final = #denovoFilter(dfinal)
    df_final = dfinal.copy()
    # df_final['INFO'] =  df_final['INFO']+';'+['denovo_filter='+val for val in df_final['denovo_filter'].to_list()]
    # TODO
    tmp_vcf = vcf_out + ".tmp"
    sample = os.path.basename(tmp_vcf).split(".")[0]

    # Create final vcf with new header
    createFinaldf(
        df_final[
            [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                sample,
            ]
        ],
        head,
        tmp_vcf,
    )

    # sort and supress dict of coverage files
    subprocess.call(
        'grep "^#" '
        + tmp_vcf
        + " > "
        + vcf_out
        + ' && grep -v "^#" '
        + tmp_vcf
        + " | sort -k1,1V -k2,2n >> "
        + vcf_out,
        shell=True,
    )
    return "ok"


df = vcfTodataframe("/STARK/output/res/MOTHER.norm.uniq.sorted.vcf")
print(df)
col = [
    "Locus_FATHER",
    "Total_Depth_FATHER",
    "Average_Depth_FATHER",
    "Depth_FATHER",
    "Base_counts_FATHER",
    "Locus_MOTHER",
    "Total_Depth_MOTHER",
    "Average_Depth_MOTHER",
    "Depth_MOTHER",
    "Base_counts_MOTHER",
    "Locus_foetus",
    "Total_Depth_foetus",
    "Average_Depth_foetus",
    "Depth_foetus",
    "Base_counts_foetus",
]
parseSpeed("/STARK/output/res/MOTHER.final.vcf.dict/COVFILEchr3", df, col, "MOTHER")
# add_info(
#    "/STARK/output/res/MOTHER.norm.uniq.sorted.vcf",
#    "/STARK/output/res/MOTHER.final.vcf",
#    "/STARK/output/res/all.gatk.cov",
# )

# origin, header = vcfTodataframe("/app/res_test/FATHER.merged.vcf", True)
# df, _, _i = parse_sample_field(origin)
##
# final = barcode_DPNI(df, "FF5")
# final.to_csv("/app/res/denovo_filter.tsv", header=True, index=False, sep="\t")
# origin["INFO"] = [
#    ";".join(item)
#    for item in list(
#        zip(
#            origin["INFO"].to_list(),
#            ["DPNI_Barcode=" + val for val in final["DPNI_Barcode"].to_list()],
#            [
#                "DPNI_Barcode_coment=" + val
#                for val in final["DPNI_Barcode_coment"].to_list()
#            ],
#            ["DPNI_Denovo=" + val for val in final["DPNI_Denovo"].to_list()],
#        )
#    )
# ]
## origin.to_csv("/app/res/test_barcode_info2.tsv", header=True, index=False, sep="\t")
# with open("/app/res/denovoinfo_filter.vcf", "w+") as f:
#    for row in header[:-1]:
#        f.write(row + "\n")
#    f.write(
#        '##INFO=<ID=DPNI_Barcode,Number=.,Type=String,Description="Family barcode from STARK DPNI module">\n'
#    )
#    f.write(
#        '##INFO=<ID=DPNI_Barcode_coment,Number=.,Type=String,Description="Comment on barcode result from STARK #DPNI module">\n'
#    )
#    f.write(
#        '##INFO=<ID=DPNI_Denovo,Number=.,Type=String,Description="denovo variants assumption">\n'
#    )
# origin.to_csv(
#    "/app/res/denovoinfo_filter.vcf", header=True, mode="a", index=False, sep="\t"
# )
