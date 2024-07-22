import os
import glob
import itertools

workdir: "/data/san/data1/users/david/mining_2023/class_IId_analysis/data/core_peptide_analyis"

model_files = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob("*.hmm")]
data_files = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob("../bakta_out/*.faa")]

def generate_combinations(models, datas):
    return expand("hmm_out/{model}_{data}.tsv", model=models, data=datas)

# Define wildcard constraints to ensure correct file matching
wildcard_constraints:
    model="|".join(model_files),
    data="|".join(data_files)

rule all:
    input:
        generate_combinations(model_files, data_files),
        "merged_hmm.tsv",
        "final_hmm_results.tsv",
        "final_proteins.faa",
        "protein_properties.tsv",
        "distinct_proteins.faa",
        "clusters/leaderless_c50id.faa",
        "clusters/leaderless_c50id.table",
        "reduced_tree/leaderless_cdhit90_aln.faa",
        "reduced_tree/leaderless_cdhit90_aln.fasta",
        "reduced_tree/leaderless_cdhit90_aln.fasta.raxml.supportFBP"

# Rule HMMs
rule hmmer:
    input:
        model = "{model}.hmm",
        data = "../bakta_out/{data}.faa"
    output:
        o1 = "hmm_out/{model}_{data}.tsv"
    threads: 16
    shell:
        """
        hmmsearch --cpu {threads} \
        --tblout {output.o1} \
        {input.model} \
        {input.data}
        """

# Rule process the hmm files to make a TSV with 
# unique IDs and <0.00001 probability
def get_tsv_files(wildcards):
    return [os.path.join(dirpath, f)
            for dirpath, dirnames, files in os.walk('hmm_out')
            for f in files if f.endswith('.tsv')]

rule concatenate_tsv:
    input:
        get_tsv_files
    output:
        'merged_hmm.tsv'
    shell:
        """
        awk 'FNR==1 && NR!=1{{next}} !/#/{{print}}' {input} > {output}
        """

# Rule to get significant hits < 0.00001
# and to keep unique IDs in column 1

rule significant_hits:
    input:
        'merged_hmm.tsv'
    output:
        'final_hmm_results.tsv'
    shell:
        """
           awk '$5 < 0.00001 && $1 != "" {print $1, $5}' {input} | awk '!seen[$0]++' > {output}
        """

# Rule get all proteins that are in column 1 of final_hmm_results.tsv
# the proteins are in bakta_out folder stored in .faa files
rule get_proteins:
    input:
        'final_hmm_results.tsv'
    output:
        'final_proteins.faa'
    conda:
        '/data/san/data1/users/david/mining_2023/class_IId_analysis/code/seqkit.yaml'
    shell:
        """
            awk '{{print $1}}' {input} > temp.txt
            while read line; do
                grep -hA 1 -w $line ../bakta_out/*.faa >> tmp2.faa
            done < temp.txt
            rm temp.txt
            cat tmp2.faa synthesized_peptides.faa mgnify_all.faa >> tmp3.faa
            seqkit rmdup -n tmp3.faa > {output}
            rm tmp2.faa tmp3.faa
        """

rule get_protein_properties:
    input:
        'final_proteins.faa'
    output:
        'protein_properties.tsv'
    conda:
        '/data/san/data1/users/david/mining_2023/class_IId_analysis/code/biopython.yml'
    shell:
        """
            python ../../code/calculate_protein_properties.py {input} {output}
        """

rule unique_sequences:
    input:
        'final_proteins.faa'
    output:
        'distinct_proteins.faa'
    conda:
        '/data/san/data1/users/david/mining_2023/class_IId_analysis/code/seqkit.yaml'
    shell:
        """
            seqkit rmdup -s {input} > {output}
        """

rule cd_hit_50:
    input:
        'final_proteins.faa'
    output:
        faa = 'clusters/leaderless_c50id.faa',
        table = 'clusters/leaderless_c50id.table'
    conda:
        '/home/david/yamls/cdhit.yml'
    shell:
        """
        cd-hit -i {input} -d 0 -n 3 -c 0.5 -o {output.faa}
        python clusters/cdhit2tbl.py clusters/leaderless_c50id.faa.clstr {output.table}
        """

rule cdhit_reduced_tree: # add an outgroup
    input:
        'final_proteins.faa'
    output:
        faa = 'reduced_tree/leaderless_cdhit90_aln.faa',
    conda:
        '/home/david/yamls/cdhit.yml'
    shell:
        """
        cd-hit -i {input} -d 0 -n 3 -c 0.9 -o {output.faa}
        cat reduced_tree/add_root.faa >> {output.faa}
        """

rule align_reduced_tree:
    input:
        faa = rules.cdhit_reduced_tree.output.faa
    output:
        aln = 'reduced_tree/leaderless_cdhit90_aln.fasta',
        fbp = "reduced_tree/leaderless_cdhit90_aln.fasta.raxml.supportFBP"
    conda:
        '/home/david/yamls/raxmlng.yml'
    shell:
        """
        muscle -in {input.faa} -out {output.aln} 
        raxml-ng --all \
            --seed 2 \
            --threads 16 \
            --msa {output.aln} \
            --model LG  \
            --bs-trees 1000 \
            --bs-metric tbe,fbp \
            --redo
        """
