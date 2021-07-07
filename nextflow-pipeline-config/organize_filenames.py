#import glob
import pandas as pd
from shutil import copyfile
from collections import defaultdict
import os

seq_dir_filenames = {"sample_name" : [], "filename0" : []}
unique = set()
source_dupes = defaultdict(list)
sample_table = pd.read_csv("sample_table.csv")
sample_files = sample_table["seq_dir"] + sample_table["fastq_filename"]
#print(sample_files)
for i, source in enumerate(sample_files):
    dest = os.path.join("manuscript_fastq", os.path.basename(source))
    if dest in unique:
        #print(source, dest)
        print("YOOOOOOOOOOOO")
        break
        dest = os.path.join("manuscript_fastq", "non-dupe-" + os.path.basename(source))
    #source_dupes[os.path.basename(source)].append(source)
    assert dest not in unique    
    unique.add(dest)
    seq_dir_filenames["sample_name"].append(f"sample_{i}")
    seq_dir_filenames["filename0"].append(f"{dest}")
    copyfile(source, dest)


#for key, value in source_dupes.items():
#    if len(value) >=2:
#        print(f"{key} has duplicates here: {value}")
#    filenames = [line.strip().split()[0] for line in open(seq_dir, "r")]
#    for f in filenames:
#        assert f not in unique
#        unique.add(f)
#    seq_dir_filenames[seq_dir[:-14] + ".tar"] = pd.Series(filenames)
seq_dir_df = pd.DataFrame(seq_dir_filenames)
#column_name_map = {i : f"filename{i}" for i in range(seq_dir_df.shape[1])}
#column_name_map["index"] = "sample_name"
#seq_dir_df.rename(column_name_map, axis=1).to_csv("sra-metadata-filenames-example.tsv", "\t")
seq_dir_df.to_csv("sra-metadata-filenames.tsv", "\t")
