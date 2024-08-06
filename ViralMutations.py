import random
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from itertools import combinations
import zipfile
import os
import zlib
from tqdm import tqdm
from math import comb
import json
import shutil

logger = logging.getLogger(__name__)
# logging.basicConfig(level=logging.INFO)

class MutationGenerator:
    # amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

    def __init__(self, args):
        self.args = args
        self.input = args.input
        self.type = args.type
        self.name = str(args.base_name)
        self.aas = 'ACDEFGHIKLMNPQRSTVWY'
        try:
            self.seq = self.read_sequence()
            logger.info("Sequence read succesful")
        except Exception as e:
            logger.error("Failed to read sequence", exc_info=True)

    def read_sequence(self):
        with open(self.input, 'r') as f:
            seq = f.readline().strip()
            seq = ''.join(line.strip() for line in f)
        logger.debug(f"Sequence loaded: {seq[:10]}...")
        return seq
    
    def load_json_file(self):
        try:
            with open(self.input, 'r') as f:
                data = json.load(f)
                logger.info("JSON file loaded succesfully")
        except Exception as e:
            logger.error("Failed to load JSON file", exc_info=True)
        return data
    
    def read_mutations(self):
        mutations = []
        with open(self.args.mutations) as f:
            for line in f:
                mutations.append(line.strip())
        return mutations

    def zip_file(self, file_counter, save_path):
        zipfile_name = os.path.join(self.args.out_dir, f'{self.args.base_name}_{file_counter}.zip')
        try:
            with zipfile.ZipFile(zipfile_name, 'w', compression=zipfile.ZIP_DEFLATED) as zipf:
                zipf.write(save_path)
            os.remove(save_path)
        except Exception as e:
            logger.error(f"Failed to zip {zipfile_name} sequence file", exc_info=True)

    def zip_all(self):
        zip_files = [f for f in os.listdir(self.args.out_dir) if f.endswith('.zip')]
        final_zip_path = os.path.join(self.args.out_dir, f"{self.args.base_name}_all_2step.zip")
        
        with zipfile.ZipFile(final_zip_path, 'w', zipfile.ZIP_DEFLATED) as final_zip:
            for zip_file in zip_files:
                zip_file_path = os.path.join(self.args.out_dir, zip_file)
                final_zip.write(zip_file_path, zip_file)
                os.remove(zip_file_path)  # Remove the original ZIP file after adding

        logger.info(f"All mutation ZIP files have been successfully combined into {final_zip_path}")

    def zip_dir(self, dir_counter, dir_path):
        zipfile_name = os.path.join(self.args.out_dir, f'{self.args.base_name}_{dir_counter}.zip')
        try:
            shutil.make_archive(base_name=zipfile_name[:-4], format='zip', root_dir=dir_path)
            shutil.rmtree(dir_path)  # Remove the directory after zipping
            logger.info(f"Zipped directory {zipfile_name}")
        except Exception as e:
            logger.error(f"Failed to zip directory {zipfile_name}", exc_info=True)

    def step_mutations(self):
        # Generate all 1-step mutations
        # if self.type == '1' or self.type == 'both':
        #     mutations = []
        #     ids = []
        #     try:
        #         for i in tqdm(range(len(self.seq)), desc='Generating 1-step sequences'):
        #             for aa in self.aas:
        #                 if aa != self.seq[i]:
        #                     mutation = self.seq[:i] + aa + self.seq[i+1:]
        #                     mutation_id = f'{self.name}-{self.seq[i]}{i+1}{aa}'
        #                     record = SeqRecord(Seq(mutation), id=mutation_id, description="")
        #                     mutations.append(record)
        #         logger.info("Successfully generated 1-step sequences")
        #     except Exception as e:
        #         logger.error("Failed to generate 1-step sequences", exc_info=True)

        #     try:
        #         save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_1step.fasta')
        #         with open(save_path, 'w') as output_handle:
        #             SeqIO.write(mutations, output_handle, "fasta")
        #         logger.info("Successfully written sequences to file")
        #     except Exception as e:
        #         logger.error("Failed to write sequences to file", exc_info=True)
        super_dir = os.path.join(self.args.out_dir, f'{self.args.base_name}_batches')
        os.makedirs(super_dir, exist_ok=True)  # Create super directory for all batches
        batch_counter = 1
        seq_counter = 0

        if self.type == '1' or self.type == 'both':
            batch_dir = os.path.join(super_dir, f'batch{batch_counter}')
            os.makedirs(batch_dir, exist_ok=True)

            try:
                for i in tqdm(range(len(self.seq)), desc='Generating 1-step sequences'):
                    for aa in self.aas:
                        if aa != self.seq[i]:
                            mutation = self.seq[:i] + aa + self.seq[i+1:]
                            mutation_id = f'{self.name}-{self.seq[i]}{i+1}{aa}'
                            record = SeqRecord(Seq(mutation), id=mutation_id, description="")
                            
                            # Writing each sequence to its own file within the batch directory
                            seq_path = os.path.join(batch_dir, f'{mutation_id}.fasta')
                            with open(seq_path, 'w') as output_handle:
                                SeqIO.write(record, output_handle, "fasta")
                            
                            seq_counter += 1
                            
                            # Check if the batch directory has 100 sequences, then move to the next batch directory
                            if seq_counter % 100 == 0:
                                batch_counter += 1
                                batch_dir = os.path.join(super_dir, f'batch{batch_counter}')
                                os.makedirs(batch_dir, exist_ok=True)

                logger.info("Successfully generated 1-step sequences")
            
            except Exception as e:
                logger.error("Failed to generate 1-step sequences", exc_info=True)
    
    # Add your 2-step mutations logic here as previously adjusted
    
    # Finally, zip the super directory containing all batches
            try:
                shutil.make_archive(base_name=os.path.join(self.args.out_dir, f'{self.args.base_name}_all_batches'), format='zip', root_dir=super_dir)
                shutil.rmtree(super_dir)  # Optional: Remove the super directory after zipping
                logger.info("Successfully zipped all batch directories")
            except Exception as e:
                logger.error("Failed to zip all batch directories", exc_info=True)
        

        if self.type == '2' or self.type == 'both':
            n = len(self.seq)
            total_seqs = int(n * 19 * (n - 1) * 19 / 2)  # Total combinations calculation
            seq_counter = 0
            dir_counter = 1
            seq_dir = os.path.join(self.args.out_dir, f'{self.args.base_name}_{dir_counter}')
            os.makedirs(seq_dir, exist_ok=True)

            try:
                with tqdm(total=total_seqs, desc='Generating 2-step sequences') as pbar:
                    for i, j in combinations(range(n), 2):
                        for aa1 in self.aas.replace(self.seq[i], ''):
                            for aa2 in self.aas.replace(self.seq[j], ''):
                                mutated_seq = list(self.seq)
                                mutated_seq[i] = aa1
                                mutated_seq[j] = aa2
                                mutated_seq = ''.join(mutated_seq)
                                record_id = f'{self.args.base_name}_{self.seq[i]}{i+1}{aa1}_{self.seq[j]}{j+1}{aa2}'
                                record = SeqRecord(Seq(mutated_seq), id=record_id, description="")
                                
                                # Writing each sequence to its own file
                                seq_path = os.path.join(seq_dir, f'{record_id}.fasta')
                                with open(seq_path, 'w') as output_handle:
                                    SeqIO.write(record, output_handle, "fasta")
                                
                                seq_counter += 1
                                pbar.update(1)
                                
                                # Check if the directory has 50 sequences, then zip and move to the next directory
                                if seq_counter % 25000 == 0:
                                    self.zip_dir(dir_counter, seq_dir)  # Zip the current directory
                                    dir_counter += 1
                                    seq_dir = os.path.join(self.args.out_dir, f'{self.args.base_name}_{dir_counter}')
                                    os.makedirs(seq_dir, exist_ok=True)
            
            except Exception as e:
                logger.error("Failed during sequence generation or file operations", exc_info=True)
            finally:
                # Handling any remaining sequences in the last directory
                if seq_counter % 10000 > 0:
                    self.zip_dir(dir_counter, seq_dir)
                logger.info(f"Successfully generated, written, and compressed {seq_counter} sequences")
        # Generate all 2-step mutations
        # if self.type == '2' or self.type == 'both':
        #     n = len(self.seq)
        #     total_seqs = (n * 19 * (n-1) * 19)/2
        #     #18915 total combs
        #     seq_counter = 0
        #     file_counter = 1
        #     save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_{file_counter}.fasta')

        #     try:
        #         with tqdm(total=total_seqs, desc=f'Generating 2-step sequences') as pbar:
        #             output_file = open(save_path, 'w')
        #             for i, j in combinations(range(n), 2):
        #                 for aa1 in self.aas.replace(self.seq[i], ''):
        #                     for aa2 in self.aas.replace(self.seq[j], ''):
        #                         mutated_seq = list(self.seq)
        #                         mutated_seq[i] = aa1
        #                         mutated_seq[j] = aa2
        #                         mutated_seq = ''.join(mutated_seq)
        #                         record = SeqRecord(Seq(mutated_seq), 
        #                                            id=f'{self.args.base_name}_{self.seq[i]}{i+1}{aa1}_{self.seq[j]}{j+1}{aa2}',
        #                                            description="")
        #                         SeqIO.write(record, output_file, 'fasta')

        #                         seq_counter += 1
        #                         pbar.update(1)
                                
        #                         if seq_counter % 1000000 == 0:
        #                             logger.info(f"The script has generated {seq_counter} sequences so far")
        #                             output_file.close()
        #                             self.zip_file(file_counter, save_path)
        #                             file_counter += 1
        #                             save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_{file_counter}.fasta')
        #                             output_file = open(save_path, 'w')

        #     except Exception as e:
        #         logger.error("Failed during sequence generation or file operations", exc_info=True)
        #     finally:
        #         output_file.close()
        #         self.zip_file(file_counter, save_path)
        #         try:
        #             self.zip_all()
        #             logger.info(f"Succesfully generated, written and compressed {seq_counter} sequences")
        #         except Exception as e:
        #             logger.error("Failed to combine all ZIP files", exc_info=True)
    

    def subsets(self, details):
        min_length, max_length, interval, repeat = details
        subset_dict = {}
        seq_length = len(self.seq)
        seq_lengths = range(min_length, max_length+1, interval)
        for length in seq_lengths:
            subset_list = []
            try:
                while len(subset_list) < repeat:
                    starting_point = random.randint(0, seq_length)
                    if starting_point - length < 0:
                        end_point = starting_point + length
                    else:
                        end_point = starting_point - length
                    start_index = min(starting_point, end_point)
                    end_index = max(starting_point, end_point)
                    subset_seq = self.seq[start_index:end_index]
                    subset_list.append(subset_seq)
                subset_dict[length] = subset_list
                logger.info(f'Generated {len(subset_list)} sequences of length {length}')
            except Exception as e:
                logger.error(f'Failed to generate {len(subset_list)} sequences of length {length}')

        save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_subset.fasta')
        try:
            with open(save_path, 'w') as output_file:
                length_list = list(subset_dict.keys())
                for length in length_list:
                    for i, sequence in enumerate(subset_dict[length]):
                        record = SeqRecord(Seq(sequence), id=f'{length}_sequence_{i+1}', description="")
                        SeqIO.write(record, output_file, 'fasta')
            logger.info("Succesfully written all subsets to file")
        except Exception as e:
            logger.error("Failed to write subsets to file")

    def specific_mutation(self, mutations):
        deletions = 0
        seq = self.seq
        for mutation in mutations:
            old_aa = mutation[0]
            new_aa = mutation[-1]
            position = int(mutation[1:-1]) - (1+deletions)
            try:
                if seq[position] == old_aa:
                    before = seq[:position]
                    after = seq[position+1:]
                    try:
                        if new_aa == '-':
                            seq = before + after 
                            logger.info(f'Position {position} has been deleted, sequence is now {len(seq)} long')
                            deletions += 1
                    except Exception as e:
                        logger.error(f"Error in deleting aa at position {position}")
            except Exception as e:
                logger.error(f'The mismatch is at {position} {seq[position]} instead of {old_aa}')
        try:
            save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_variant.fasta')
            with open(save_path, 'w') as f:
                record = SeqRecord(Seq(seq), id=f'{self.args.base_name}_variant', description="")
                SeqIO.write(record, f, "fasta")
            logger.info("Succesfully implemented mutations and written to file")
        except Exception as e:
            logger.error("Succesfully implemented mutations but failed to write to file")

    def write_from_json(self):
        try:
            data = self.load_json_file(self)
            save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_json.fasta')
            for key in data.keys():
                with open(save_path, 'w') as f:
                    record = SeqRecord(Seq(data[key]), id=key, description="")
                    SeqIO.write(record, f, "fasta")
            logger.info("Succesfully written sequences to file")
        except Exception as e:
            logger.error("Failed to write sequences to file", exc_info=True)
    
    def alternative_mutations(self, mutations):
        alt_mutations = []
        try:
            for mut in mutations:
                old_change = mut[-1:]
                pos = int(mut[1:4])-1
                new_change = mut[0]
                for aa in self.aas:
                    if aa != old_change:
                        mutation = self.seq[:pos] + aa + self.seq[pos+1:]
                        mutation_id = f'{self.name}-{new_change}{pos}{aa}'
                        record = SeqRecord(Seq(mutation), id=mutation_id, description="")
                        alt_mutations.append(record)
            logger.info("Successfully generated alternative mutations")
        except Exception as e:
            logger.error("Failed to generate alternative mutations", exc_info=True)
        try:
            save_path = os.path.join(self.args.out_dir, f'{self.args.base_name}_alternatives.fasta')
            with open(save_path, 'w') as output_handle:
                SeqIO.write(alt_mutations, output_handle, "fasta")
            logger.info("Successfully written sequences to file")
        except Exception as e:
            logger.error("Failed to write sequences to file", exc_info=True)