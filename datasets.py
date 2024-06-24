from dataclasses import dataclass
from pathlib import Path
from enum import IntEnum, StrEnum

import numpy as np
import pandas as pd


SNP_SV_DIR = Path("data", "snp-sv")
TRAIT_DIR = Path("data", "traits")
GENE_DIR = Path("data", "genes")

class ViewType(IntEnum):
    SNP = 0
    SV = 1
    Trait = 2
    Region = 3


@dataclass
class StructuralVariant:
    sv_id: str
    start: int
    end: int
    center: int
    length: int
    variant: str
    chromosome: int

    @classmethod
    def parse(cls, sv_id: str) -> "StructuralVariant":
        try: 
            chr_part, rest = sv_id.split(":")
            start, rest = rest.split("-")
            end, variant = rest.split("_")
            chr_num = int(chr_part[3:])
            length = int(end) - int(start)
            center = int(end) - int((int(length) / 2))
            return cls(sv_id=sv_id, start=start, end=end, center=center, length=length, variant=variant, chromosome=chr_num)
        except ValueError:
            chr_part, pos, type, subtype = sv_id.split("_")
            chr_num = int(chr_part[3:])
            return cls(sv_id=sv_id, start=pos, end=int(pos)+1, center=pos, length=1, variant=f"{type}_{subtype}", chromosome=chr_num)
    

@dataclass(init=False)
class Data():
    snp_sv_df: pd.DataFrame
    sv_df: pd.DataFrame
    trait_df: pd.DataFrame
    snps: list[str]
    svs: list[StructuralVariant]
    traits: list[str]

    def __init__(self, trait: str, p_value: float = 0.00005):
        snp_sv_df = pd.read_csv(Path(SNP_SV_DIR, "chr1.allQTLs.NEWSET.JOIN_size.txt"), delimiter="\t")
        snp_sv_df = snp_sv_df[["#CHROM", "POS", "REF", "ALT", "T_STAT", "P", "SV_ID"]].rename(columns={"#CHROM": "CHR", "SV_ID": "sv_id"})
        
        trait = "alzheimer"
        trait_df = pd.read_csv(Path(TRAIT_DIR, trait, "chr1_Alzheimer_million_hg38.txt"), delimiter="\t")
        trait_df = trait_df[["CHR", "POS", "P", "RSID"]]

        snps = [snp for snp in snp_sv_df["POS"].drop_duplicates()]
        sv_ids = [sv_id for sv_id in snp_sv_df["sv_id"].drop_duplicates()]
        svs = [StructuralVariant.parse(sv_id) for sv_id in sv_ids]
        traits = [trait for trait in trait_df["RSID"].drop_duplicates()]

        # bonferroni_threshold_snp_sv_df = p_value / len(snp_sv_df)
        # snp_sv_df = snp_sv_df.loc[snp_sv_df["P"] < bonferroni_threshold_snp_sv_df]
        snp_sv_df["P"] = snp_sv_df["P"].astype(float)
        snp_sv_df["-log10P"] = snp_sv_df["P"].apply(lambda x: -np.log10(float(x)))

        sv_df = pd.DataFrame([sv.__dict__ for sv in svs])

        # bonferroni_threshold_trait_df = p_value / len(trait_df)
        # trait_df = trait_df.loc[trait_df["P"] < bonferroni_threshold_trait_df]
        trait_df["P"] = trait_df["P"].astype(float)
        trait_df["-log10P"] = trait_df["P"].apply(lambda x: -np.log10(float(x)))

        self.snp_sv_df = snp_sv_df
        self.sv_df = sv_df
        self.trait_df = trait_df
        self.snps = snps
        self.svs = svs
        self.traits = traits


data = Data("alzheimer")
