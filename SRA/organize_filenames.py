#import glob
import pandas as pd
from shutil import copyfile
from collections import defaultdict
import os
import sys

#sra_meta = {"sample_name" : [], "filename0" : []}
sra_meta = defaultdict(list)
bio_attr = defaultdict(list)
unique = set()
source_dupes = defaultdict(list)
sample_table = pd.read_csv(sys.argv[1])
#sample_files = sample_table["seq_dir"] + sample_table["fastq_filename"]
#print(sample_files)
# for i, source in enumerate(sample_files):
for i, (row, data) in enumerate(sample_table.iterrows()):
    #source = data["seq_dir"] + data["fastq_filename"]
    #dest = os.path.join("manuscript_fastq", os.path.basename(source))
    #if dest in unique:
        #print(source, dest)
        #print("YOOOOOOOOOOOO")
        #break
        #dest = os.path.join("manuscript_fastq", "non-dupe-" + os.path.basename(source))
    #source_dupes[os.path.basename(source)].append(source)
    #assert dest not in unique    
    #unique.add(dest)
    bio_attr["sample_name"].append(f"sample_{i}")
    sra_meta["sample_name"].append(f"sample_{i}")
    sra_meta["filename0"].append(os.path.basename(data["fastq_filename"]))
    lib_id = "-".join(data["experiment"].split("/"))
    sra_meta["library id"].append(data["library_batch"]+"-"+lib_id)
    bio_attr["collection date"].append(data["experiment"])

bio_attr["isolate"]="not applicable"
bio_attr["host"]="Escherichia coli"
bio_attr["lab host"]="Escherichia coli"
bio_attr["geographic location"]="USA:Seattle"
bio_attr["isolation source"] = "T7 phage library"

sra_meta["Organism:"] = "Escherichia phage T7"
sra_meta["Instrument:"] = "Illumina MiSeq"
sra_meta["Strategy:"] = "CLONEEND"
sra_meta["Source:"] = "GENOMIC"
sra_meta["Selection:"] = "PCR"
sra_meta["Layout:"] = "SINGLE"
    #copyfile(source, dest)


#for key, value in source_dupes.items():
#    if len(value) >=2:
#        print(f"{key} has duplicates here: {value}")
#    filenames = [line.strip().split()[0] for line in open(seq_dir, "r")]
#    for f in filenames:
#        assert f not in unique
#        unique.add(f)
#    sra_meta[seq_dir[:-14] + ".tar"] = pd.Series(filenames)
pd.DataFrame(sra_meta).set_index("sample_name").to_csv("sra_meta.tsv", "\t")
pd.DataFrame(bio_attr).set_index("sample_name").to_csv("bio_attr.tsv", "\t")
#column_name_map = {i : f"filename{i}" for i in range(seq_dir_df.shape[1])}

#column_name_map["index"] = "sample_name"
#seq_dir_df.rename(column_name_map, axis=1).to_csv("sra-metadata-filenames-example.tsv", "\t")
