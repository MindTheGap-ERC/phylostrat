NUMBER_OF_REPLICAS=1

BASE_COLLECTION="irods://nluu11p/home/research-mindthegap/phylostrat"
BASE_COLLECTION="./outputs"
REPLICA_IDS = range(NUMBER_OF_REPLICAS)
Id = [1]


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
        
ID = ['1']
NCHAR = ['30']

rule create_replicas:
    input:
        char=expand("data/sim/char_mat_cont_rho0_nchar{nchar}_{id}.nex", id = ID, nchar = NCHAR),
        fossils=expand("data/sim/fossils_cont_rho0_{id}.tsv", id = ID, nchar = NCHAR)
    output:
        #storage(BASE_COLLECTION + "/simulation/mcmc_trace_{run_id}.log")
        numeric=expand("output/res_cont_{id}_{nchar}.log", id = ID, nchar = NCHAR),
        trees=expand("output/trees_base_cont_{id}_{nchar}.log", id = ID, nchar = NCHAR)
    shell:
        "bash code/inference_test.sh {input.char} {input.fossils} {output.numeric} {output.trees}"

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