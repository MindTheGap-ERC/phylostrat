set -euo pipefail
set -x

INPUT_MORPH=$1
INPUT_TAXA=$2
INPUT_MOL=$3
OUTPUT_NUMERIC=$4
OUTPUT_TREES=$5

RB_COMMAND="morph <- readDiscreteCharacterData(\"${INPUT_MORPH}\"); taxa <- readTaxonData(\"${INPUT_TAXA}\"); mol <- readDiscreteCharacterData(\"${INPUT_MOL}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\"; source(\"code/under_prior/fbd_inference_skyline_A_te_nd.rev\")" 

echo $RB_COMMAND | "rb"
echo $?
echo "te nd skyline inference done"
