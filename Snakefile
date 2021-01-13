"""
Binning MAGs for John
Hey John
You probably already know what my whole process is but I guess this is just my
notes and commentary on the small refining steps I like to do. 
"""

# execution params
max_threads = config.get('max_threads', 9)

###########
# Set Up Variables

min_contig_length = int(config.get('min_contig_length', 1500))
use_finishm = config.get('use_finishm', True)
mag_dir = "MAGs_finished" if use_finishm else "MAGs"

# Set up input files
#  This should be a dict of sample -> {"fastq": reads_fastq, "contigs": fasta_file}
fasta_template = "assembly/{sample}/contigs.all.fasta"
fastq_template = "assembly/{sample}/{assembly}.clean.fastq"
samples, assemblies = glob_wildcards(fastq_template)
logger.debug("Found {} samples".format(len(samples)))
assembly_files = {s:{"contigs": fasta_template.format(sample=s),
                     "fastq": fastq_template.format(sample=s, assembly=a),
                     "bam": f"binning/{s}/coverm/contigs.all.fasta.{a}.clean.fastq.bam"}
                  for s, a in zip(samples, assemblies)}

# binning tools
method_params_list = \
    [f'maxbin.{min_contig_length}.{markers}' for markers in ['40', '107']] + \
    [f'metabat2.{min_contig_length}', 'concoct.default'] + \
    [f'metabat.{min_contig_length}.{mode}' \
     for mode in ['specific', 'veryspecific', 'sensitive', \
                  'verysensitive', 'superspecific']] 
logger.debug("METHODS: " + repr(method_params_list))

###########
# Target
#  the first rule defines the default output
localrules: output

#output_files = expand("binning/{sample}/das_tool_DASTool_scaffolds2bin.txt", \
#                      sample=samples)
output_files = [f'{mag_dir}/dRep'] 
logger.debug("OUTPUT FILES: " + repr(output_files))

rule output:
    input: output_files

rule up_to_checkpoints:
    input:
        expand("binning/{sample}/das_tool_DASTool_summary.txt",
               sample=samples)

rule through_checkpoints:
    input:
        expand("binning/{sample}/merged_bins",
               sample=samples)

###########
# Coverm
#        So starting from mapping I like to use a program called coverM 
#        (its great cause it uses minimap2 and is fast at calculating a 
#        overage file for metabat).
 
rule coverm_make:
    """
    Map: generate bam files mapping reads to contigs for each sample
    """
    input:
        fastq=lambda w: assembly_files[w.sample]['fastq'],
        contigs=lambda w: assembly_files[w.sample]['contigs']
    output:
        bam="binning/{sample}/coverm/{contigs}.fasta.{reads}.fastq.bam"
    threads: max_threads
    params:
        dir="binning/{sample}/coverm"
    shell: 
       "coverm make --interleaved {input.fastq} --reference {input.contigs} \
                    --threads {threads} --output-directory {params.dir} \
            > {output.bam}.log 2>&1"

rule coverm_contig:
    """ Generate counts """
    input: lambda w: assembly_files[w.sample]['bam']
    output: "binning/{sample}/coverm.{sample}.metabat.tsv"
    threads: max_threads
    shell:
        "coverm contig -t {threads} -b {input} -m metabat \
            > {output} \
            2> {output}.log"

rule coverm_gawk:
    """ reformat coverm counts """
    input: "binning/{sample}/coverm.{sample}.metabat.tsv"
    output: "binning/{sample}/coverm.{sample}.abund.tsv"
    shell:
        """gawk '{{print $1"\t"$3}}' {input} > {output}"""

###########
# Binning
#    metabat, metabat2, maxbin, concoct

rule metabat_sensitive:
    """ metabat is run a few times. The first time we save the TNF file for reuse """
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs'],
        counts="binning/{sample}/coverm.{sample}.metabat.tsv"
    output:
        bins=directory("binning/{sample}" \
                       "/metabat.{min_contig}.sensitive.bins"),
        tnf="binning/{sample}/metabat.{min_contig}.tnf",
        dist="binning/{sample}/metabat.{min_contig}.dist",
    threads: max_threads
    shell: """
        rm -rf {output.bins}
        metabat1 -t {threads} -i {input.contigs} -a {input.counts} \
                     --minContig {wildcards.min_contig} -v -B 20 --keep \
                     --saveTNF {output.tnf} --saveDistance {output.dist} \
                     --sensitive \
                     -o {output.bins}/{wildcards.sample}.bin \
            > {output.bins}.log 2>&1
        """
    
rule metabat_other:
    """ metabat is run a few more times. the mode comes from das_tool inputs """
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs'],
        counts="binning/{sample}/coverm.{sample}.metabat.tsv",
        tnf="binning/{sample}/metabat.{min_contig}.tnf",
        dist="binning/{sample}/metabat.{min_contig}.dist",
    output:
        bins=directory("binning/{sample}" \
                       "/metabat.{min_contig}.{mode}.bins"),
    wildcard_constraints:
        mode="(specific|veryspecific|superspecific|verysensitive)"
    threads: max_threads
    shell: """
        rm -rf {output.bins}
        metabat1 -t {threads} -i {input.contigs} -a {input.counts} \
                     --minContig {wildcards.min_contig} -v -B 20 --keep \
                     --saveTNF {input.tnf} --saveDistance {input.dist} \
                     --{wildcards.mode} \
                     -o {output.bins}/{wildcards.sample}.bin \
            > {output.bins}.log 2>&1
        """
   
rule metabat2:
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs'],
        counts="binning/{sample}/coverm.{sample}.metabat.tsv"
    output:
        bins=directory("binning/{sample}" \
                       "/metabat2.{min_contig}.bins"),
    threads: max_threads
    shell: """
        rm -rf {output.bins}
        metabat2 -t {threads} -i {input.contigs} -a {input.counts} \
                     --minContig {wildcards.min_contig} \
                     -o {output.bins}/bin \
            > {output.bins}.log 2>&1
        """

rule maxbin:
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs'],
        counts="binning/{sample}/coverm.{sample}.abund.tsv"
    output:
        bins=directory("binning/{sample}" \
                       "/maxbin.{min_contig}.{markers}.bins"),
    params:
        marker_set=lambda w: "" if w.markers == '107' \
                             else f"-markerset {w.markers}"
    threads: max_threads
    shell: """
        rm -rf {output.bins}
        mkdir -p {output.bins}
        run_MaxBin.pl -min_contig_length {wildcards.min_contig} \
                          -thread {threads} -contig {input.contigs} \
                          -out {output.bins}/bin \
                          -abund {input.counts} {params.marker_set} \
            > {output.bins}.log 2>&1
       """

# concoct binning happens in a few steps
rule concoct_cut_up:
    """ Concoct binning: cup up fasta into chunks """
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs']
    output:
        fasta="binning/{sample}/concoct_files/contigs_10K.fa",
        bed="binning/{sample}/concoct_files/contigs_10K.bed"
    shell: "cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last \
                            -b {output.bed} > {output.fasta} \
                2> {output.fasta}.log"

rule index_bam:
    input: "{file_root}.bam"
    output: "{file_root}.bam.bai"
    threads: max_threads
    shell: "samtools index {input} -@ {threads} > {output}.log 2>&1"

rule concoct_coverage_table:
    input:
        bed="binning/{sample}/concoct_files/contigs_10K.bed",
        bai_files=lambda w: assembly_files[w.sample]['bam'] + ".bai"
    output: "binning/{sample}/concoct_files/coverage_table.tsv"
    params:
        bam_files=lambda w: assembly_files[w.sample]['bam']
    shell:
        "concoct_coverage_table.py {input.bed} {params.bam_files} \
          > {output} 2> {output}.log"

rule concoct:
    input:
        fasta="binning/{sample}/concoct_files/contigs_10K.fa",
        coverage="binning/{sample}/concoct_files/coverage_table.tsv"
    output:
        "binning/{sample}/concoct_files/concoct_output_clustering_gt1000.csv"
    params: 
        out_prefix="binning/{sample}/concoct_files/concoct_output"
    threads: max_threads
    shell:
        "concoct -t {threads} --composition_file {input.fasta} \
                 --coverage_file {input.coverage} -b {params.out_prefix} \
            > {output}.log 2>&1 "

rule concoct_merge:
    input: rules.concoct.output
    output: "binning/{sample}/concoct_files/concoct.clustering_merged.csv"
    shell:
        "merge_cutup_clustering.py {input} > {output} 2> {output}.log"

rule concoct_extract_fasta:
    input:
        contigs=lambda w: assembly_files[w.sample]['contigs'],
        clusters=rules.concoct_merge.output
    output: directory("binning/{sample}/concoct.default.bins")
    shell:
        """rm -rf {output}
           mkdir -p {output}
           extract_fasta_bins.py {input.contigs} {input.clusters} \
                               --output_path {output} \
             > {output}.log 2>&1
        """

##########
# Merge Bins wit DAS_tool
rule make_bin_table:
    "turn folder of MAG fasta files into table mapping file name to contig name"
    input: "binning/{sample}/{method}.{params}.bins"
    output: "binning/{sample}/{method}.{params}.scaffolds2bin.tsv"
    wildcard_constraints:
        method="(metabat2?|maxbin|concoct)"
    threads: 3
    shell: """( grep -H "^>" {input}/*.f*a || true ) \
               | sed 's/:>/\t/g' \
               | gawk '{{print $2"\t"$1}}' \
               > {output}"""

rule das_tool:
    input:
        files=expand("binning/{{sample}}/{method_params}.scaffolds2bin.tsv",
                     method_params=method_params_list),
        contigs=lambda w: assembly_files[w.sample]['contigs'],
    output:
        bins="binning/{sample}/das_tool_DASTool_scaffolds2bin.txt",
        summary="binning/{sample}/das_tool_DASTool_summary.txt"
    threads: max_threads
    params:
        n=len(method_params_list),
        out_pref="binning/{sample}/das_tool",
        files=lambda w: ",".join( \
            (f"{w.sample}/binning/{method_params}.scaffolds2bin.tsv" \
             for method_params in method_params_list)),
        labels=",".join(method_params_list),
    shell:
        """
        # clean up from any previous runs
        rm -rf {params.out_pref}*

        # Remove empty files from the inputs
        L={params.n}
        N=0
        INPUT_FILES=({input.files})
        LABELS=({method_params_list})
        INPUT_FILES_FILE={output.bins}.infiles
        LABEL_FILE={output.bins}.labels

        while [ "$N" -lt "$L" ]; do
            if [ -s ${{INPUT_FILES[$N]}} ]; then
                if [ -s ${{INPUT_FILES_FILE}} ]; then
                    # Insert a comma if there is already data in the file
                    echo -n , >> $INPUT_FILES_FILE
                    echo -n , >> $LABEL_FILE
                fi
                echo -n ${{INPUT_FILES[$N]}} >> $INPUT_FILES_FILE
                echo -n ${{LABELS[$N]}} >> $LABEL_FILE
            fi
            let N=N+1
        done
        # Pull file contests into variables
        INPUT_FILES=$(cat $INPUT_FILES_FILE)
        LABELS=$(cat $LABEL_FILE)

        # now run das_tool
        DAS_Tool -i ${{INPUT_FILES}} -l ${{LABELS}} -c {input.contigs} \
                  -o {params.out_pref} -t {threads} --search_engine diamond \
            > {output.summary}.log 2>&1
        """

"""
Part 2 of the workflow reassembles the DAS_Tool bins, so we'll need a
checkpoint rule to extract the bin to fasta files and then the rest will go
from there


DAS_Tool can extract the bins for us, but it doesn't like the files being in a
subdirectory. 
"""

checkpoint das_tool_bins:
    """
     Convert das_tool output to fasta files. Das_Tool can do this,
     but it doesn't work if the input fasta files are in subdirs.
     
     Also, remove _sub from some fasta file names. I'm not sure why
     das_tool is adding that.
     """
    input: 
        "binning/{sample}/das_tool_DASTool_scaffolds2bin.txt",
    output: directory("binning/{sample}/merged_bins")
    shell:
        """
        rm -rf {output}
        mkdir -p {output}
        cut -f 2 {input} | uniq | while read FASTA; do
          FASTA=${{FASTA%%_sub}}
          METHOD=$(basename $(dirname $FASTA))
          METHOD=${{METHOD%%.bins}}
          OUTF={output}/$METHOD.$(basename $FASTA)
          OUTL={output}/$METHOD.$(basename $FASTA).list
          grep "$FASTA" {input} | cut -f 1 > $OUTL
          seqtk subseq $FASTA $OUTL > $OUTF 2> $OUTF.log
        done
        """

def get_sample_bins(sample):
    checkpoint_dir = checkpoints.das_tool_bins.get(sample=sample).output
    bins, = glob_wildcards(f"binning/{sample}/merged_bins/{{bin_id}}.fasta")
    return bins

rule extract_reads:
    """ 
     two steps (1: get read names from hits 2: pull reads by name) so we get
     pairs even if only one matched
    """
    input:
        contigs="binning/{sample}/merged_bins/{bin_id}.fasta.list",
        bai_files=lambda w: assembly_files[w.sample]['bam'] + ".bai"
    output: "reassembly/{sample}/{bin_id}.reads.fastq"
    params:
        bam_files=lambda w: assembly_files[w.sample]['bam']
    threads: 2
    shell: """
        rm -f {output} {output}.list {output}.bed

        # convert contig list to BED
        cat {input.contigs} | while read C; do
            echo -e "$C\t0\t9999999" ;
        done > {output}.bed

        # use bash piped input to get read list and fastq from bam file
        for BAM_FILE in {params.bam_files}; do
            seqtk subseq \
                <(samtools fastq -n $BAM_FILE) \
                <(samtools view -L {output}.bed $BAM_FILE | cut -f 1) \
                >> {output}
        done
        """

rule reassemble_spades:
    input:
        reads="reassembly/{sample}/{bin_id}.reads.fastq",
        prev_contigs="binning/{sample}/merged_bins/{bin_id}.fasta"
    output:
        "reassembly/{sample}/{bin_id}.spades/scaffolds.fasta"
    params:
        out_dir="reassembly/{sample}/{bin_id}.spades"
    threads: max_threads
    conda: "conda.spades.yaml"
    shell:
        """
        rm -rf {params.out_dir}
        spades.py -t {threads} --trusted-contigs {input.prev_contigs} \
          --12 {input.reads} -o {params.out_dir} \
          > {params.out_dir}.log 2>&1
        """

rule min_contig_length:
    input:
        raw="reassembly/{sample}/{bin_id}.spades/scaffolds.fasta"
    output:
        filtered="reassembly/{sample}/{bin_id}.filtered.fasta"
    run:
        from Bio import SeqIO
        with open(output.filtered, 'wt') as out_handle:
            for contig in SeqIO.parse(input.raw, 'fasta'):
                if len(contig) >= int(min_contig_length):
                    out_handle.write(contig.format('fasta'))

rule finishm:
    """
    Running 11 of these in parallel killed our system,
    so it may be wise to use the made up finishm resource
    to limit the number that run.
    """
    input:
        fasta="reassembly/{sample}/{bin_id}.filtered.fasta",
        fastq=lambda w: assembly_files[w.sample]['fastq'],
    output:
        fasta="reassembly/{sample}/finished/{bin_id}.filtered.fasta.scaffolds.fasta"
    conda: "conda.finishm.yaml"
    resources:
        finishm=1   
    params:
        out_dir="reassembly/{sample}/finished"
    shell: """
        mkdir -p {params.out_dir}
        finishm roundup --genomes {input.fasta} \
                --interleaved-fastq {input.fastq} \
                --output-directory {params.out_dir} \
            > {output.fasta}.log 2>&1
        """

rule link_for_checkm:
    """
    checkm needs all the genomes in one folder, so let's use symlinks to collect 
    """
    input: rules.finishm.output if use_finishm else rules.reassemble_spades.output
    output: "{mag_dir}/{sample}_{bin_id}.fasta"
    shell: "ln -s ../{input} {output}"

rule checkm:
    input:
        fastas=lambda w: \
                   [f"{mag_dir}/{sample}_{bin_id}.fasta" \
                    for sample in samples \
                    for bin_id in get_sample_bins(sample) \
                   ],
    output: "{mag_dir}/checkm.tsv",
    threads: max_threads
    shell: """
        checkm lineage_wf ./{mag_dir} -x fasta ./{mag_dir}/checkm --tab_table -t \
            {threads} -f {output} > {output}.log 2>&1
        """

rule checkm_o2:
    input: rules.checkm.output
    output: "{mag_dir}/checkm.o2.tsv"
    shell:
        "checkm qa {mag_dir}/checkm/lineage.ms {mag_dir}/checkm -o2 --tab_table \
            -f {output} > {output}.log 2>&1"

rule genome_info:
    """ reformat the CheckM output, also force second checkm to run """
    input:
        checkm=rules.checkm.output,
        o2=rules.checkm_o2.output
    output: "{mag_dir}/genome_info.csv"
    shell:
        """
        gawk -F"\t" '{{print $1".fasta,"$12","$13}}' {input.checkm} \
         | sed 's/Bin Id.fasta,Completeness,Contamination/genome,completeness,contamination/g' \
         > {output} 2> {output}.log
        """

rule dereplicate:
    input:
        fastas=lambda w: \
                   [f"{mag_dir}/{sample}_{bin_id}.fasta" \
                    for sample in samples \
                    for bin_id in get_sample_bins(sample) \
                   ],
        info=rules.genome_info.output
    output: directory("{mag_dir}/dRep")
    conda: "conda.drep.yaml"
    threads: max_threads
    shell: """
        rm -rf {output}
        dRep dereplicate {output} -g {input.fastas} -p {threads} \
          --genomeInfo {input.info} -comp 70 -con 10 -nc 0.6 -sa 0.97 \
          > {output}.log 2>&1
        """
