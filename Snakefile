BASE_COLLECTION="irods://nluu11p/home/research-mindthegap/phylostrat"

SIM_LOCATION = "data/sim"
MCMC_LOCATION = "output/"
NO_OF_REPLICAS = 50
#replica ID
ID = [str(i) for i in range(1, NO_OF_REPLICAS + 1)]
NCHAR = ["30", "300", "1000"] # ["30", "300", "1000"]
CASE =  ["cont", "inc_A", "inc_B"] #["cont"] #
RHO = ["0", "1"]
NRUNS_NUM = 2 # number of independent chains run per analysis. must be larger than 2
NRUNS = str(NRUNS_NUM)
NRUN = [str(i) for i in range(1, NRUNS_NUM+1)]

# rule all:
#     input:
#         numeric_cont=expand("output/res_const_{case}_{id}_{nchar}_run_{run}.log", id = ID, nchar = NCHAR, case = CASE, run = NRUN),
#         trees_cont=expand("output/trees_const_{case}_{id}_{nchar}_run_{run}.log", id = ID, nchar = NCHAR, case = CASE, run = NRUN),
#         numeric_skyline_A=expand("output/res_skyline_A_{id}_{nchar}_run_{run}.log", id = ID, nchar = NCHAR, run = NRUN),
#         trees_skyline_A=expand("output/trees_skyline_A_{id}_{nchar}_run_{run}.log", id = ID, nchar = NCHAR, run = NRUN)
#         #storage(expand(BASE_COLLECTION+"/simulation/mcmc_trace_{run_id}.log", run_id=REPLICA_IDS))

storage yoda:
    provider="irods"

rule all:
    input:
        storage.yoda(expand(BASE_COLLECTION + "/simulation_data/tree_complete_{id}.nex", id = ID)),
        storage.yoda(expand(BASE_COLLECTION + "/simulation_data/mol_data_{id}.nex", id = ID)),
        storage.yoda(expand(BASE_COLLECTION + "/simulation_data/tree_all_fossils_{id}.nex", id = ID)),
        storage.yoda(expand(BASE_COLLECTION + "/simulation_data/fossils_{case}_rho{rho}_{id}.tsv", id = ID, rho = RHO, case = CASE)),
        storage.yoda(expand(BASE_COLLECTION + "/simulation_data/char_mat_{case}_rho{rho}_nchar{nchar}_{id}.nex", case = CASE, rho = RHO, nchar = NCHAR, id = ID))


rule copy_data_to_yoda:
    output:
        full_tree = storage.yoda(BASE_COLLECTION + "/simulation_data/tree_complete_{id}.nex"),
        mol = storage.yoda(BASE_COLLECTION + "/simulation_data/mol_data_{id}.nex"),
        fossil_tree = storage.yoda(BASE_COLLECTION + "/simulation_data/tree_all_fossils_{id}.nex")       
    input:
        full_tree = "data/sim/tree_complete_{id}.nex",
        mol = "data/sim/mol_data_{id}.nex",
        fossil_tree = "data/sim/tree_all_fossils_{id}.nex"        
    shell:
        """
        mv {input.full_tree} {output.full_tree}
        mv {input.mol} {output.mol}
        mv {input.fossil_tree} {output.fossil_tree}
        """

rule copy_fossils_to_yoda:
    output:
        fossils = storage.yoda(BASE_COLLECTION + "/simulation_data/fossils_{case}_rho{rho}_{id}.tsv")
    input:
        fossils = "data/sim/fossils_{case}_rho{rho}_{id}.tsv"
    shell:
        "mv {input.fossils} {output.fossils}"

rule copy_char_to_yoda:
    output:
        char_mat = storage.yoda(BASE_COLLECTION + "/simulation_data/char_mat_{case}_rho{rho}_nchar{nchar}_{id}.nex")
    input:
        char_mat = "data/sim/char_mat_{case}_rho{rho}_nchar{nchar}_{id}.nex"
    shell:
        "mv {input.char_mat} {output.char_mat}"

rule simulate_data_locally:
    output:
        expand("data/sim/tree_complete_{id}.nex", id = ID),
        expand("data/sim/fossils_{case}_rho{rho}_{id}.tsv", case = CASE, id = ID, rho = RHO),
        expand("data/sim/tree_all_fossils_{id}.nex", id = ID),
        expand("data/sim/mol_data_{id}.nex", id = ID),
        expand("data/sim/char_mat_{case}_rho{rho}_nchar{nchar}_{id}.nex", case = CASE, rho = RHO, nchar = NCHAR, id = ID)
    shell:
        "Rscript code/sim.R"




# tree_comp = storage(expand("irods://nluu11p/home/research-mindthegap/phylostrat/simulation_data/tree_complete_{id}.nex", id = ID))
# rule fbd_analysis_constant_rates:
#     params:
#         nruns = NRUNS,
#         numeric="output/res_const_{case}_{id}_{nchar}.log",
#         trees="output/trees_const_{case}_{id}_{nchar}.log"
#     input:
#         char="data/sim/char_mat_{case}_rho0_nchar{nchar}_{id}.nex",
#         fossils="data/sim/fossils_{case}_rho0_{id}.tsv"
#     output:
#         #storage(BASE_COLLECTION + "/simulation/mcmc_trace_{run_id}.log")
#         numeric=expand("output/res_const_{{case}}_{{id}}_{{nchar}}_run_{run}.log", run = NRUN),
#         trees=expand("output/trees_const_{{case}}_{{id}}_{{nchar}}_run_{run}.log", run = NRUN)
#     shell:
#         "bash code/inference_fbd_const_rate.sh {input.char} {input.fossils} {params.numeric} {params.trees} {params.nruns}"

# rule fbd_analysis_skyline_A:
#     params:
#         nruns = NRUNS,
#         numeric="output/res_skyline_A_{id}_{nchar}.log",
#         trees="output/trees_skyline_A_{id}_{nchar}.log"
#     input:
#         char = "data/sim/char_mat_inc_A_rho0_nchar{nchar}_{id}.nex",
#         fossils="data/sim/fossils_inc_A_rho0_{id}.tsv"
#     output:
#         numeric=expand("output/res_skyline_A_{{id}}_{{nchar}}_run_{run}.log", run = NRUN),
#         trees=expand("output/trees_skyline_A_{{id}}_{{nchar}}_run_{run}.log", run = NRUN)
#     shell:
#         "bash code/inference_fbd_skyline_A.sh {input.char} {input.fossils} {params.numeric} {params.trees} {params.nruns}"

# rule run_model:
#     # input:
#         # param_file = storage(BASE_COLLECTION+ "params.json")
#     output:
#         csv = storage(BASE_COLLECTION + "/results/output_{run_id}.csv"),
#         toml = storage(BASE_COLLECTION + "/results/output_{run_id}.toml"),
#         h5 = storage(BASE_COLLECTION + "/results/output_{run_id}.h5")    
#     shell:
#         """
#         julia --project=. run_dissolution.jl param_dissolution.json {wildcards.run_id} output_{wildcards.run_id}.csv output_{wildcards.run_id}.toml output_{wildcards.run_id}.h5
                
#         cp output_{wildcards.run_id}.csv {output.csv}
#         cp output_{wildcards.run_id}.toml {output.toml}
#         cp output_{wildcards.run_id}.h5 {output.h5}
        
#         rm output_{wildcards.run_id}.csv output_{wildcards.run_id}.toml output_{wildcards.run_id}.h5
#         """

# rule generate_plot:
#     input:
#         result_file = storage(BASE_COLLECTION+"/results/output_{run_id}.csv")
#     output:
#         storage(BASE_COLLECTION+"/plots/plot_{run_id}.png")
#     shell:
#         """
#         julia plot_results.jl {input} {output}
#         """