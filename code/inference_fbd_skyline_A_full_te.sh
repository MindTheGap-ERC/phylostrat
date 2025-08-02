INPUT_MORPH=$1
INPUT_TAXA=$2
INPUT_MOL=$3
OUTPUT_NUMERIC=$4
OUTPUT_TREES=$5

RB_COMMAND="morph <- readDiscreteCharacterData(\"${INPUT_MORPH}\"); taxa <- readTaxonData(\"${INPUT_TAXA}\"); mol <- readDiscreteCharacterData(\"${INPUT_MOL}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\"; source(\"code/fbd_inference_skyline_A_full_te.rev\")" 

echo $RB_COMMAND | "rb"