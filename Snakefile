import json

NUMBER_OF_REPLICAS=50

BASE_COLLECTION="irods://nluu11p/home/research-mindthegap/phylostrat"
REPLICA_IDS = range(NUMBER_OF_REPLICAS)


rule all:
    input:
        storage(expand(BASE_COLLECTION+"/simulation/mcmc_trace_{run_id}.log", run_id=REPLICA_IDS))
        # storage(expand(BASE_COLLECTION+"/plots/plot_{run_id}.png", run_id=ID))

# rule create_params:
#     output:
#         param_result = storage(BASE_COLLECTION+ "params.json")
#     shell:
#         """julia --project=. param_input.jl {output.param_result} """
        

rule create_replicas:
    output:
        storage(BASE_COLLECTION + "/simulation/mcmc_trace_{run_id}.log")

    shell:
        "bash code/inference_base.sh {wildcards.run_id} 1000 {output}"

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