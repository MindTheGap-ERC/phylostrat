set -euo pipefail
INPUT_CHAR=$1
INPUT_FOSSILS=$2
OUTPUT_NUMERIC=$3
OUTPUT_TREES=$4
INPUT_MOL=$5
rb_command="morph_name <- \"${INPUT_CHAR}\"; taxa_name <- \"${INPUT_FOSSILS}\";  outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\";  mol_name <- \"${INPUT_MOL}\"; source(\"code/fbd_inference_informative_sl.rev\")"
echo $rb_command | "rb"