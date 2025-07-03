NUMBER_OF_REPLICAS=1

BASE_COLLECTION="irods://nluu11p/home/research-mindthegap/phylostrat"
BASE_COLLECTION="./outputs"


#rule all:
#    input:
#        "output/res_base_cont_30_1.log"
        #storage(expand(BASE_COLLECTION+"/simulation/mcmc_trace_{run_id}.log", run_id=REPLICA_IDS))
        # storage(expand(BASE_COLLECTION+"/plots/plot_{run_id}.png", run_id=ID))

# rule create_params:
#     output:
#         param_result = storage(BASE_COLLECTION+ "params.json")
#     shell:
#         """julia --project=. param_input.jl {output.param_result} """
        
ID =   ["1"]#[str(i) for i in range(1, 2)]
NCHAR = ["30"] # ["30", "300", "1000"]
CASE = ["cont"] # ["cont", "inc_A", "inc_B"]

rule fbd_analysis_constant_rates:
    input:
        char=expand("./data/sim/char_mat_{case}_rho0_nchar{nchar}_{id}.nex", id = ID, nchar = NCHAR, case = CASE),
        fossils=expand("./data/sim/fossils_{case}_rho0_{id}.tsv", id = ID, case = CASE)
    output:
        #storage(BASE_COLLECTION + "/simulation/mcmc_trace_{run_id}.log")
        numeric=expand("./output/res_const_{case}_{id}_{nchar}.log", id = ID, nchar = NCHAR, case = CASE),
        trees=expand("./output/trees_const_{case}_{id}_{nchar}.log", id = ID, nchar = NCHAR, case = CASE)
    shell:
        "bash code/inference_test.sh {input.char} {input.fossils} {output.numeric} {output.trees}"

rule fbd_analysis_skyline_A:
    input:
        char = expand("data/sim/char_mat_inc_A_rho0_nchar{nchar}_{id}.nex", id = ID, nchar = NCHAR),
        fossils=expand("data/sim/fossils_inc_A_rho0_{id}.tsv", id = ID)
    output:
        numeric=expand("output/res_skyline_A_{id}_{nchar}.log", id = ID, nchar = NCHAR),
        trees=expand("output/trees_skyline_A_{id}_{nchar}.log", id = ID, nchar = NCHAR)
    shell:
        "bash code/inference_fbd_skyline_A.sh {input.char} {input.fossils} {output.numeric} {output.trees}"

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