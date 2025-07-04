INPUT_CHAR=$1
INPUT_FOSSILS=$2
OUTPUT_NUMERIC=$3
OUTPUT_TREES=$4
NRUNS=$5
rb_command="morph <- readDiscreteCharacterData(\"${INPUT_CHAR}\"); taxa <- readTaxonData(\"${INPUT_FOSSILS}\"); outfile <- \"${OUTPUT_NUMERIC}\"; treefile_name <- \"${OUTPUT_TREES}\"; nruns <- \"${NRUNS}\"; source(\"code/skyline_A_inference.rev\")"
echo $rb_command | "rb"