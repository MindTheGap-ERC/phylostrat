INPUT_CHAR=$1
INPUT_FOSSILS=$2
OUTPUT_NUMERIC=$3
OUTPUT_TREES=$4
#NCHAR=$2
#OUTPUT_TREE_FILE=$3
#INS="cont"
#echo $ID
#OUTPUT_NUMERIC_FILE="output/res_base_${INS}_${ID}_${NCHAR}.log"
#echo $OUTPUT_NUMERIC_FILE
#echo $INS
#rb_command="morph <- readDiscreteCharacterData(\"data/sim/char_mat_${INS}_rho0_nchar${NCHAR}_${ID}.nex\"); taxa <- readTaxonData(\"data/sim/fossils_${INS}_rho0_${ID}.tsv\"); outfile <- \"${OUTPUT_NUMERIC_FILE}\"; treefile_name <- \"${OUTPUT_TREE_FILE}\"; source(\"code/fbd_inference.rev\")"
rb_command="morph <- readDiscreteCharacterData(\"${INPUT_CHAR}\"); taxa <- readTaxonData(\"${INPUT_FOSSILS}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile_name <- \"${OUTPUT_TREES}\"; source(\"code/fbd_inference.rev\")"
echo $rb_command | "rb.exe"