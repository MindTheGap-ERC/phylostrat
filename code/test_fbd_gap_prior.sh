ID=1
NCHAR=30
CASE=miller
bash code/inference_fbd_gap_prior.sh data/sim/fbd_strat/char_mat_${ID}_nchar${NCHAR}_${CASE}.nex data/sim/fbd_strat/fossils_${ID}_nchar${NCHAR}_${CASE}.csv data/sim/fbd_strat/tree_${ID}_nchar${NCHAR}_${CASE}.nex output/fbd_strat/num_test_${ID}_${NCHAR}_${CASE} output/fbd_strat/tree_test_${ID}_${NCHAR}_${CASE} data/sim/fbd_strat/mol_dat_${ID}_nchar${NCHAR}_${CASE}.nex