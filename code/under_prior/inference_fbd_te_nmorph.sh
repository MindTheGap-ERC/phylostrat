set -euo pipefail
set -x

# fbd inference - total evidence constant rates
INPUT_MORPH=$1
INPUT_TAXA=$2
INPUT_MOL=$3
OUTPUT_NUMERIC=$4
OUTPUT_TREES=$5

RB_COMMAND="morph <- readDiscreteCharacterData(\"${INPUT_MORPH}\"); taxa <- readTaxonData(\"${INPUT_TAXA}\"); mol <- readDiscreteCharacterData(\"${INPUT_MOL}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\"; source(\"code/under_prior/fbd_inference_te_nmorph.rev\")" 

echo $RB_COMMAND | "rb"
echo $?
echo "te nmorph inference done"
