echo "bovine STAB"
cd 1fq9_bovine_stab/charmm-gui-8035848484/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_bovine_stab/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_bovine_stab/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_bovine_stab/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_bovine_stab/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_bovine_stab/compare_total_rep5.png
cd ../../../
echo "Goldfish"
cd 1fq9_goldfish/charmm-gui-8035857277/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_goldfish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_goldfish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_goldfish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_goldfish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_goldfish/compare_total_rep5.png
cd ../../../
echo "human"
cd 1fq9_human/charmm-gui-8035859451/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human/compare_total_rep5.png
cd ../../../
echo "human STAB"
cd 1fq9_human_stab/charmm-gui-8035865058/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human_stab/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human_stab/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human_stab/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human_stab/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_human_stab/compare_total_rep5.png
cd ../../../
echo "Japanese Ricefish"
cd 1fq9_japanese_ricefish/charmm-gui-8042199431/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_japanese_ricefish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_japanese_ricefish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_japanese_ricefish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_japanese_ricefish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_japanese_ricefish/compare_total_rep5.png
cd ../../../
echo "Zebrafish"
cd 1fq9_zebrafish/charmm-gui-8035871929/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_zebrafish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_zebrafish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_zebrafish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_zebrafish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1fq9_zebrafish/compare_total_rep5.png
cd ../../../
