#bash code/inference_fbd_const_rate_te.sh data/sim/char_mat_cont_rho1_nchar30_1.nex data/sim/fossils_cont_rho1_1.tsv data/sim/mol_data_1.nex output/test1.log output/test_tree1.log
#bash code/inference_fbd_skyline_A.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/res_skyline_A_1_30.log output/trees_skyline_A_1_30.log
#bash code/inference_fbd_skyline_A_te.sh data/sim/char_mat_cont_rho1_nchar30_1.nex data/sim/fossils_cont_rho1_1.tsv data/sim/mol_data_1.nex output/test1.log output/test_tree1.log

#bash code/inference_fbd_fo.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/res_skyline_A_1_30.log output/trees_skyline_A_1_30.log

#bash code/inference_fbd_skyline_A_fo.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/res_skyline_A_1_30.log output/trees_skyline_A_1_30.log

bash code/inference_fbd_fo.sh data/sim/char_mat_cont_rho0_nchar30_1.nex data/sim/fossils_cont_rho0_1.tsv output/num_test_fo_cont.log output/tree_test_fo_cont.log &
bash code/inference_fbd_fo.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/num_test_fo_inc_A.log output/tree_test_fo_inc_A.log &
bash code/inference_fbd_fo.sh data/sim/char_mat_inc_B_rho0_nchar30_1.nex data/sim/fossils_inc_B_rho0_1.tsv output/num_test_fo_inc_B.log output/tree_test_fo_inc_B.log &
bash code/inference_fbd_skyline_A_fo.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/num_test_fo_skyline_A.log output/tree_test_fo_skyline_A.log &
bash code/inference_fbd_skyline_A_full_fo.sh data/sim/char_mat_inc_A_rho0_nchar30_1.nex data/sim/fossils_inc_A_rho0_1.tsv output/num_test_fo_skyline_A_full.log output/tree_test_fo_skyline_A_full.log &
wait