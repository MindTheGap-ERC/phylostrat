ID=$1
NCHAR=$2
OUTPUT_TREE_FILE=$3
INS="cont"
OUTPUT_NUMERIC_FILE="output/res_base_${INS}_${NCHAR}_${ID}.log":
rb_command="morph <- readDiscreteCharacterData(\"data/sim/char_mat_${INS}_rho0_nchar${NCHAR}_${ID}.nex\"); taxa <- readTaxonData(\"data/sim/fossils_${INS}_rho0_${ID}.tsv\"); outfile <- \"${OUTPUT_NUMERIC_FILE}\"; treefile_name <- \"${OUTPUT_TREE_FILE}\"; source(\"code/fbd_inference.rev\")"

echo $rb_command | "rb"