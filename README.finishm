== Installation ==
Simply create the conda environment specified in `snakemake.yaml`:

    $ mamba env create -p ../mag.env -f ./snakemake.yaml

=== Databases ===
Checkm2: See note below about the CHECKM2 database. It's easy to install.

GTDBtk: You'll need to download the whole thing and unpack it. It's frickin huge. (TODO: instructions)

== Running ==
The workflow requires a set of read fastq files and contig fasta files from a
collection of metagenomic short read assemblies. These should be named
consistently so that they can be found from templates. The template parsing
assumes one set of contigs per assembly and possibly multiple sets of reads per
assembly.

I recommend running outside of the repo dir which requires specifying the
makefile path with `-s`.

An example command:

    $ snakemake -s /path/to/mags_from_separate_assemblies/Snakefile \
        --config \
        fasta_template=/path/to/assemblies/{assembly}/contigs.fasta \
        fastq_tempalte=/path/to/assemblies/{assembly}/reads.{read_set}.fastq \
        gtdbtk_dir=/path/to/GTDBtk/data \
        -j 20 --use-conda --conda-frontend mamba

=== Additional parameters ===

Additional configuration parameters:

 * max_threads (9): no single task can request more than this number of threads (except concoct)
 * concoct_threads (12): threads for the concoct tasks
 * checkm2_db_path (None): location of the checkm2 database file (optional. See note below)
 * min_contig_length (1500): ignore contigs shorter than this
 * use_finishm (True): set to False to skip the finishm step (sometimes the conda installation fails and it's not strictly necessary)
 * skip_samples (): pass a semi-colon-separated list of assemblies to ignore if your templates pull in too many.

== NOTES ==

=== CheckM2 ===
Checkm2 is installed by conda, but you need a database. Fortunately, this is
easy. Once you have created the conda environment from `snakemake.yaml`, run
the following command in the environment (to install in your home folder):

    $ checkm2 database --download

If you want to put the DBN somewhere else, use the `--path` option:

    $ checkm2 database --download --path /put/it/somehwere/else

If you already downloaded it with a previous checkm2 installation, you can just
set the path with:

    $ checkm2 database --setdblocation /existing/db/path/uniref100.KO.1.dmnd

Be aware that for the `--path` argument, you are giving it a directory to create
and put the DB file in. For the `--setdblocation` argument, you are passing in
the full name of the diamond dartabase file, not just the folder.

Alternatively, you can pass the db location at runtime with the snakemake
configuration parameter `checkm2_db_path`. EG:

    $ snakemake --config checkm2_db_path=/existing/db/path/uniref100.KO.1.dmnd ...

=== FINISHM ===
The finishm step of this workflow is fragile.

I have created my own conda package to install finishmto try to solve this with limited success. On some systems, I have had to go in and manually recompile the velvet binaries installed by ruby into the conda env.

    $ cd ${CONDA_ENV}/lib/ruby/gems/2.5.0/gems/finishm-0.0.9/ext/src
    $ make -j2 velvetg velveth
    
If this doesn't work or finishim is too slow, you can skip it in the workflow by setting the use_finishm config value to False:

    $ snakemake --config use_finishm=False ...
    
