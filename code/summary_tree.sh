set -euo pipefail
TREEFILE=$1
rb_command="treefile <- \"${TREEFILE}\"; source(\"code/summary_tree.rev\")"
echo $rb_command | "rb"