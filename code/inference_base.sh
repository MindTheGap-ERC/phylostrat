ID=$1
NCHAR=$2
INS="cont"

rb_command="morph <- readDiscreteCharacterData(\"data/sim/char_mat_${INS}_rho0_nchar${NCHAR}_${ID}.nex\"); taxa <- readTaxonData(\"data/sim/fossils_${INS}_rho0_${ID}.tsv\"); outfile <- \"output/res_base_${INS}_${NCHAR}_${ID}.log\"; treefile_name <- \"output/trees_base_${INS}_${NCHAR}_${ID}.log\"; source(\"code/fbd_inference.rev\")"

echo $rb_command | "rb.exe"