# fossil only fbd inference - skyline A model
set -euo pipefail
INPUT_CHAR=$1
INPUT_FOSSILS=$2
OUTPUT_NUMERIC=$3
OUTPUT_TREES=$4
rb_command="morph <- readDiscreteCharacterData(\"${INPUT_CHAR}\"); taxa <- readTaxonData(\"${INPUT_FOSSILS}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile_name <- \"${OUTPUT_TREES}\"; source(\"code/fbd_inference_skyline_A_fo.rev\")"
echo $rb_command | "rb"