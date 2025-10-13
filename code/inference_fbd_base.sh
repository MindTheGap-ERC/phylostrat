set -euo pipefail
INPUT_CHAR=$1
INPUT_FOSSILS=$2
INPUT_TREE=$3
OUTPUT_NUMERIC=$4
OUTPUT_TREES=$5
INPUT_MOL=$6
rb_command="morph_name <- \"${INPUT_CHAR}\"; taxa_name <- \"${INPUT_FOSSILS}\"; tree_name <- \"${INPUT_TREE}\";  outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\";  mol_name <- \"${INPUT_MOL}\"; source(\"code/fbd_inference_base.rev\")"
echo $rb_command | "rb"