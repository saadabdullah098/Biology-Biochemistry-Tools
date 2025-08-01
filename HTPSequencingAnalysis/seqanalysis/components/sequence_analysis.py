import os
import re
from glob import glob
from typing import List, Tuple, Dict
import hashlib
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from .benchling_connection import BenchlingConnection
from seqanalysis import constants


class SequenceAnalysis:
    """
    The Sequence analysis tool compares vendor sequencing output to reference files.
    The vendor output is adjusted for orientation and origin to compare with the reference
    Exact length and sequence matches are determined and reported

    parameters
    ----------
    - date
    - sample directory
    - reference directory
    
    returns
    -------
    output csv file with Pass/Fail calls
    """

    def __init__(self,
                 date: str = None,
                 sample_directory: str = None):
        self.date = date
        self.sample_dir = sample_directory

    def compare_sequences(self):
        """
            Compares the reference and sequenced sequence
        """
        comparison = self.shift_sequenced_plasmid()
        comparison = self.set_parent_hash(comparison)
        comparison = self.set_sequenced_hash(comparison)

        seq_pass = list()
        length_pass = list()
        rows = comparison.index.tolist()
        for val in rows:
            if comparison.at[val, "parent_hash"] == comparison.at[val, "cloned_hash"]:
                seq_pass.append("Pass")
            else:
                seq_pass.append("Fail")

            if comparison.at[val, "parent_length"] == comparison.at[val, "shifted_seq_length"]:
                length_pass.append("Pass")
            else:
                length_pass.append("Fail")

        comparison["length_match"] = length_pass
        comparison["sequence_match"] = seq_pass
        comparison.rename(columns={"parent_length": "expected_size", "shifted_seq_length": "actual_size"}, inplace=True)
        print(f"final_df is {comparison}")
        comparison = comparison[["reference_plasmid", "cloned_plasmid", "well", "sequencing_dimer",
                                 "length_match", "sequence_match", "expected_size", "actual_size"]]
        comparison.to_csv(f"~/Desktop/{self.date}_Plasmid_Sequence_Comparison.csv", index=False)

        return comparison

    @staticmethod
    def set_sequenced_hash(comparison) -> pd.DataFrame:
        """
            Sets sha256 hash on the shifter upper case sequenced plasmid sequence
        """
        seq_hash = list()
        rows = comparison.index.tolist()
        for val in rows:
            seq = comparison.at[val, "shifted_sequence"].upper().encode()
            hash_out = hashlib.sha256(seq)
            seq_hash.append(hash_out.hexdigest())

        comparison["cloned_hash"] = seq_hash

        return comparison

    @staticmethod
    def set_parent_hash(comparison) -> pd.DataFrame:
        """
            Sets sha256 hash on upper case parent plasmid sequence
        """
        seq_hash = list()
        rows = comparison.index.tolist()
        for val in rows:
            seq = comparison.at[val, "parent_sequence"].upper().encode()
            hash_out = hashlib.sha256(seq)
            seq_hash.append(hash_out.hexdigest())

        comparison["parent_hash"] = seq_hash

        return comparison

    def shift_sequenced_plasmid(self) -> pd.DataFrame:
        """
            Compare sequenced plasmid to reference. restore origin and orientation for sha256 comparison
        """
        comparison = self.join_seq_to_ref()

        seq_shift = list()
        shift_length = list()
        rows = comparison.index.tolist()
        for val in rows:
            if "dimer" in comparison.at[val, "sequencing_dimer"]:
                seq_cut = int(len(comparison.at[val, "raw_sequence"]) / 2)
                print(seq_cut)
                start = comparison.at[val, "seq_delim"]
                target = comparison.at[val, "raw_sequence"][:seq_cut]
                target_rc = str(Seq(target).reverse_complement())

                match = re.search(start, target)
                rc_match = re.search(start, target_rc)

                if match is not None:
                    first = target[match.span()[1]:]
                    second = target[: match.span()[0]]
                    new_seq = start + first + second
                    seq_shift.append(new_seq)
                    shift_length.append(len(new_seq))

                elif rc_match is not None:
                    first = target_rc[rc_match.span()[1]:]
                    second = target_rc[: rc_match.span()[0]]
                    new_seq = start + first + second
                    seq_shift.append(new_seq)
                    shift_length.append(len(new_seq))
                else:
                    seq_shift.append("fail")
                    shift_length.append(0)
            else:
                start = comparison.at[val, "seq_delim"]
                target = comparison.at[val, "raw_sequence"]
                target_rc = str(Seq(comparison.at[val, "raw_sequence"]).reverse_complement())

                match = re.search(start, target)
                rc_match = re.search(start, target_rc)

                if match is not None:
                    first = target[match.span()[1]:]
                    second = target[: match.span()[0]]
                    new_seq = start + first + second
                    seq_shift.append(new_seq)
                    shift_length.append(len(new_seq))

                elif rc_match is not None:
                    first = target_rc[rc_match.span()[1]:]
                    second = target_rc[: rc_match.span()[0]]
                    new_seq = start + first + second
                    seq_shift.append(new_seq)
                    shift_length.append(len(new_seq))
                else:
                    seq_shift.append("fail")
                    shift_length.append(0)

        comparison["shifted_sequence"] = seq_shift
        comparison["shifted_seq_length"] = shift_length

        return comparison

    def join_seq_to_ref(self) -> pd.DataFrame:
        """
            Join reference data to sequenced data for comparison
        """
        seq_files = self.import_files()
        seqs = self.process_sequence_files(seq_files)

        seqs["base_plasmid_name"] = seqs["cloned_plasmid"].apply(str).str.split("-").str[0]
        unique_base_plasmids = seqs["base_plasmid_name"].unique().tolist()

        refs = self.get_reference_DNA(unique_base_plasmids)
        refs["base_plasmid_name"] = refs["reference_plasmid"].apply(str).str.split("_").str[0]

        comparison_base = seqs.merge(refs, on="base_plasmid_name", how="inner")
        comparison = self.dimer_check(comparison_base)
        comparison["seq_delim"] = comparison["parent_sequence"].str[:300]

        return comparison

    @staticmethod
    def dimer_check(comparison) -> pd.DataFrame:
        """
            Checks sequencing output length against designed lengths to identify monomers, dimers, and indels
        """
        dimer_check = list()
        rows = comparison.index.tolist()
        for val in rows:
            ratio = (comparison.at[val, "cloned_length"] / comparison.at[val, "parent_length"])
            if ratio == 1.0:
                dimer_check.append("monomer")
            elif ratio == 2.0:
                dimer_check.append("dimer")
            else:
                dimer_check.append("indel")
        comparison["sequencing_dimer"] = dimer_check

        return comparison

    @staticmethod
    def get_reference_DNA(base_plasmids) -> pd.DataFrame:
        """
            Fetch the DNA sequence from Benchling for each ID in a given list
        """
        con = BenchlingConnection("production")
        benchling, connection = con.benchling_connection()

        all_data = []

        # Process each DNA sequence ID in the list
        for dna_id in base_plasmids:
            # Fetch the DNA sequence for the current ID
            dna_sequence = benchling.dna_sequences.list(name=[dna_id])

            # Iterate through pages and sequences
            for page in dna_sequence:
                for sequence in page:
                    # Calculate the length of the bases
                    base_length = len(sequence.bases)
                    # Append the name, bases, and base length to the data list
                    all_data.append({
                        "reference_plasmid": sequence.name,
                        "parent_sequence": sequence.bases.upper(),
                        "parent_length": base_length
                    })

        refs = pd.DataFrame(all_data)

        return refs

    def process_sequence_files(self,
                               seq_files) -> pd.DataFrame:
        """
            Extract sequence, name, and sequence length from sequences genbank files
        """
        col_dict = self.column_dict()
        sequence = list()
        plas = list()
        length = list()
        wells = list()
        for file in seq_files:
            rec = SeqIO.parse(file, "genbank")
            for record in rec:
                sequence.append(str(record.seq).upper())
                plas.append(str(record.name))
                length.append(len(str(record.seq)))
            well = col_dict[int(file.split("_")[-2])]
            wells.append(well)
        seqs = pd.DataFrame()
        seqs["cloned_plasmid"] = plas
        seqs["raw_sequence"] = sequence
        seqs["cloned_length"] = length
        seqs["well"] = wells


        return seqs


    def import_files(self) -> Tuple[List, List]:
        """
            Imports a list of files from the sample directory
        """
        seq_path = f"{self.sample_dir}"
        seq_string = os.path.expanduser(seq_path)
        seq_files = [file for path, subdir, files in os.walk(seq_string) for file in glob(os.path.join(path, "*.gb*"))]

        return seq_files

    def column_dict(self) -> Dict:
        """
            Translates numbers to 96 well plate well format
        """
        col_dict = constants.well_col_dict

        return col_dict

