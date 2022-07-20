import pandas as pd
from tqdm import tqdm
import os


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


origin, header = vcfTodataframe("/app/res_test/FATHER.merged.vcf", True)
df, _, _i = parse_sample_field(origin)
#
final = barcode_DPNI(df, "FF5")
final.to_csv("/app/res/denovo_filter.tsv", header=True, index=False, sep="\t")
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
