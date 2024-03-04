#!/bin/sh

tar xfz aln_data.tar.gz
tar xfz cnt_data.tar.gz
gunzip -k pmam2023_05_output_gc211004.gz

stat='max'

./fig1_percid_gaps3.R --pub --Bthresh=0.9 --Cgapl=5 --pdf fig1mam_yb_i_${stat} --stat g_${stat} -t mam5_yst1_bct1.times -f aln_data/human_vs_mam_canon.summ2i,aln_data/mouse_vs_rat_canon.summ2i,aln_data/ecoli_vs_bct_canon.summ2i,aln_data/yeast_vs_ystc_canon.summ2i
mv fig1mam_yb_i_${stat}_id_gaps_pub.pdf figs/

## figure 2 -- pipeline

## figure 3 -- fig3_multi_tree.pdf

## figures 4,5

./fig5_percid_gaps6.R --pub --supp --pdf fig5y_i_${stat} --stat g_${stat} -t mam5_yst1_bct1.times -Y -f aln_data/prop_HUMMUS_v_all.summ2i_m,aln_data/canon_HUMMUS_v_all.summ2i_m,aln_data/conf_HUMMUS_v_all.summ2i_m,aln_data/yeast_vs_ystc_canon.summ2i_m,aln_data/ecoli_vs_bct_canon.summ2i_m
mv fig5y_i_${stat}_id_gapsF3_main_pub.pdf figs/fig5_canon_prop_pub.pdf
mv fig5y_i_${stat}_id_gapsF3_supp_pub.pdf figs/fig4_conf_pub.pdf

## figure 6
./fig6_plot_changes.R --pub --pdf fig6 -f history_new_changes_byspecies.tsv
mv fig6_cumm_prop_pub.pdf figs/

## figure 7

./fig7_clade_stats.R --pub --pdf fig7 -f cnt_data/qfomam_2022_05_clade_stats2.tab,cnt_data/mam_2023_01_clade_stats2.tab,cnt_data/mam_2023_02_clade_stats2.tab,cnt_data/mam_2023_03_clade_stats2.tab,cnt_data/mam_2023_04_clade_stats2.tab,cnt_data/mam_2023_05_clade_stats2.tab
mv fig7_clade_size_pub.pdf figs/

## suppl fig 1
./sfig1_cost_hist.R --pub -f qfomam_output_confirm240206,qfomam_output_changes240206.one
mv supp_fig1_cost_pub.pdf figs/

## suppl fig 3
sfig3_plen_diff.R --pub --pdf supp_fig3 -f pmam2023_05_output_gc211004
mv supp_fig3_ldelta_pub.pdf figs/

## suppl fig 4
./fig1_percid_gaps3.R --pub --Bthresh=0.7 --Cgapl=10 --pdf supp_fig4_mamvrtpln_i_${stat}_id70-10 --stat g_${stat} -t mam1_pln2_vrt1_yst1_bct1.times -f aln_data/arath_vs_pln2_canon.summ2i,aln_data/human_vs_mouse_10090.summ2i,aln_data/human_vs_xentr_8364_VT40.summ2i,aln_data/ecoli_vs_bct_canon.summ2i,aln_data/yeast_vs_ystc_canon.summ2i
mv supp_fig4_mamvrtpln_i_max_id70-10_id_gaps_pub.pdf figs/
