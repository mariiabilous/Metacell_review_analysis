from pathlib import Path

gamma_values = [75]
samples_list = list(range(1, 51))
datasets_list = [3,7,10,12,15,20,25,30,35,40,50]
tools = ["SuperCell", "SEACells", "MetaCell"]

annotation_labels = ["embryo_id"]
seurat_inputs = ["sc_BPCells", "sc"]

rule all:
    input:
        'results/data_full/sc_data.rds',
        'results/data_full/metadata.csv',
        'results/data_full/gene.csv',
         expand('results/data_split/sc_data/sc_data_noBPcells_{sample}.rds', sample=samples_list),
         expand('results/data_split/sc_data/sc_data_BPcells_{sample}.rds', sample=samples_list),
         expand('results/MC_anno_{annotation}/{tool}/gamma{gamma}/dataset{dataset}/merged_metacells.rds', dataset=datasets_list, gamma=gamma_values, annotation=annotation_labels, tool = tools),
         expand('results/MC_anno_{annotation}/benchmarks_time/{input_type}_{dataset}datasets_downstream.txt', dataset=datasets_list, annotation=annotation_labels, input_type = seurat_inputs),
         expand('results/MC_anno_{annotation}/benchmarks_time/MC_gamma{gamma}_{tool}_{dataset}datasets_downstream.txt', dataset=datasets_list, gamma=gamma_values, annotation=annotation_labels, tool = tools)

rule dataset_download:
    '''
    Download HLCA dataset.
    '''
    output:
        sc_data = 'results/data_full/sc_data.rds',
        cell_annotation = 'results/data_full/metadata.csv',
        genes_annotation = 'results/data_full/gene.csv'
    container:
        "docker://agabriel/matk:gpu"
    params:
        download_link_data = config['data']['data_link_mouse'],
        download_link_metadata = config['data']['annotation_link_mouse'],
        download_link_genes = config['data']['genes_link_mouse']
    shell:
        """
        wget -O results/data_full/sc_data.rds {params.download_link_data};
        wget -O results/data_full/metadata.csv {params.download_link_metadata};
        wget -O results/data_full/gene.csv {params.download_link_genes}
        """

rule dataset_splitting:
    '''
    Generate datasets with different sizes.
    '''
    input:
        sc_data = 'results/data_full/sc_data.rds',
        cell_annotation = 'results/data_full/metadata.csv',
        genes_annotation = 'results/data_full/gene.csv'
    output:
        sc_data = 'results/data_split/sc_data/sc_data_noBPcells_{sample}.rds',
        sc_data_bpcells = 'results/data_split/sc_data/sc_data_BPcells_{sample}.rds'
    container:
        "docker://agabriel/matk:gpu"
    log:
        'logs/dataset_splitting/dataset_splitting_{sample}.txt'
    shell:
        """
        Rscript scripts/dataset_splitting.r {input.sc_data} {input.genes_annotation} {input.cell_annotation} {wildcards.sample} results/data_split/sc_data/ > {log} 2>&1
        """

rule run_MATK_per_sample:
    '''
    Build metacells for each data data subset and using each metacell construction tool.
    '''
    input:
        sc_data = expand('results/data_split/sc_data/sc_data_noBPcells_{sample}.rds', sample=samples_list)
    output:
        mc_data = 'results/MC_anno_{annotation}/{tool}/gamma{gamma}/dataset{dataset}/merged_metacells.rds',
        MetaCell_benchmarks = 'results/MC_anno_{annotation}/benchmarks_time/gamma{gamma}_{dataset}_{tool}.txt'
    container:
        "docker://agabriel/matk:gpu"
    benchmark:
        'results/MC_anno_{annotation}/benchmarks/gamma{gamma}_{dataset}_{tool}.txt'
    log:
        'results/MC_anno_{annotation}/logs/run_MATK_per_sample/gamma{gamma}_{dataset}_{tool}.txt'
    threads: 15
    resources:
        mem_mb = config["run_MATK_per_sample"]["mem_mb"],
        time = config["run_MATK_per_sample"]["time"],
        partition = lambda wildcards: config[wildcards.tool]["partition"],
        gpus = lambda wildcards: config[wildcards.tool]["gpus"]
    shell:
        """
        mkdir -p results/MC_anno_{wildcards.annotation}/benchmarks_time;
        /usr/bin/time -v -o results/MC_anno_{wildcards.annotation}/benchmarks_time/gamma{wildcards.gamma}_{wildcards.dataset}_{wildcards.tool}.txt \
            Rscript scripts/run_MATK_per_sample.r -t {wildcards.tool} -d {wildcards.dataset} -g {wildcards.gamma} -n {threads} -o results/MC_anno_{wildcards.annotation}/{wildcards.tool}/gamma{wildcards.gamma}/dataset{wildcards.dataset}/ > {log} 2>&1
        """

rule run_MATK_per_sample_downstream:
    '''
    Build metacells for each data data subset and using each metacell construction tool + downstream.
    '''
    input:
        sc_data = expand('results/data_split/sc_data/sc_data_noBPcells_{sample}.rds', sample=samples_list)
    output:
        mc_integration = 'results/MC_anno_{annotation}/{tool}_downstream/gamma{gamma}/dataset{dataset}/datasets_integrated_umap.pdf',
        MetaCell_downstream_benchmarks = 'results/MC_anno_{annotation}/benchmarks_time/MC_gamma{gamma}_{tool}_{dataset}datasets_downstream.txt'
    container:
        "docker://agabriel/matk:gpu"
    benchmark:
        'results/MC_anno_{annotation}/benchmarks/gamma{gamma}_{dataset}_{tool}_downstream.txt'
    log:
        'results/MC_anno_{annotation}/logs/run_MATK_per_sample_downstream/gamma{gamma}_{dataset}_{tool}_downstream.txt'
    threads: 15
    resources:
        mem_mb = config["run_MATK_per_sample_downstream"]["mem_mb"],
        time = config["run_MATK_per_sample_downstream"]["time"],
        partition = lambda wildcards: config[wildcards.tool]["partition"],
        gpus = lambda wildcards: config[wildcards.tool]["gpus"]
    wildcard_constraints:
        tool="|".join(tools)
    shell:
        """
        mkdir -p results/MC_anno_{wildcards.annotation}/benchmarks_time;
        /usr/bin/time -v -o results/MC_anno_{wildcards.annotation}/benchmarks_time/MC_gamma{wildcards.gamma}_{wildcards.tool}_{wildcards.dataset}datasets_downstream.txt \
            Rscript scripts/MC_seurat_downstream_analysis.r -t {wildcards.tool} -d {wildcards.dataset} -g {wildcards.gamma} -n {threads} -o results/MC_anno_{wildcards.annotation}/{wildcards.tool}_downstream/gamma{wildcards.gamma}/dataset{wildcards.dataset}/ > {log} 2>&1
        """

rule run_downstream_analysis:
    '''
    Generate datasets with different sizes.
    '''
    input:
        anndata_file_sc1 = expand('results/data_split/sc_data/sc_data_BPcells_{sample}.rds', sample=samples_list),
        anndata_file_sc2 = expand('results/data_split/sc_data/sc_data_noBPcells_{sample}.rds', sample=samples_list)

    output:
        pdf_umap = 'results/MC_anno_{annotation}/downstream_sc/{input_type}_{dataset}datasets_integrated_umap.pdf',
        benchmarking_res = 'results/MC_anno_{annotation}/benchmarks_time/{input_type}_{dataset}datasets_downstream.txt'
    container:
        "docker://agabriel/matk:gpu"
    resources:
        mem_mb = config['run_downstream_analysis']['mem_mb'],
        time = config['run_downstream_analysis']['time']
    benchmark:
        'results/MC_anno_{annotation}/benchmarks/{input_type}_{dataset}datasets_downstream.txt'
    log:
        'results/MC_anno_{annotation}/logs/run_downstream_analysis/{input_type}_{dataset}datasets_downstream.txt'
    wildcard_constraints:
        input_type="|".join(seurat_inputs)
    shell:
        """
        /usr/bin/time -v -o results/MC_anno_{wildcards.annotation}/benchmarks_time/{wildcards.input_type}_{wildcards.dataset}datasets_downstream.txt \
            Rscript scripts/seurat_downstream_analysis.r {wildcards.dataset} {wildcards.input_type} results/data_split/sc_data/ results/MC_anno_{wildcards.annotation}/downstream_sc/ > {log} 2>&1
        """
