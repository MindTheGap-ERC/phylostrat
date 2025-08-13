INPUT_MORPH=$1
INPUT_TAXA=$2
INPUT_MOL=$3
OUTPUT_NUMERIC=$4
OUTPUT_TREES=$5
START_TREE=$6

RB_COMMAND="morph_name <- \"${INPUT_MORPH}\" ; taxa_name <- \"${INPUT_TAXA}\"; mol_name <- \"${INPUT_MOL}\"; outfile <- \"${OUTPUT_NUMERIC}\"; treefile <- \"${OUTPUT_TREES}\"; start_tree <- \"${START_TREE}\"; source(\"code/te/fbd_inference_skyline_A_te.rev\")" 


echo $RB_COMMAND | "rb"