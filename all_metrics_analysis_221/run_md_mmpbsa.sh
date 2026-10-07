echo "bovine STAB"
cd 1e0o_bovine_stab/charmm-gui-8219507732/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_bovine_stab/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_bovine_stab/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_bovine_stab/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_bovine_stab/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_bovine_stab/compare_total_rep5.png
cd ../../../
echo "Goldfish"
cd 1e0o_goldfish/charmm-gui-8219546278/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_goldfish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_goldfish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_goldfish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_goldfish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_goldfish/compare_total_rep5.png
cd ../../../
echo "human"
cd 1e0o_human/charmm-gui-8224696110/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human/compare_total_rep5.png
cd ../../../
echo "human STAB"
cd 1e0o_human_stab/charmm-gui-8225124362/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human_stab/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human_stab/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human_stab/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human_stab/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_human_stab/compare_total_rep5.png
cd ../../../
echo "Japanese Ricefish"
cd 1e0o_japanese_ricefish/charmm-gui-8225128173/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_japanese_ricefish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_japanese_ricefish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_japanese_ricefish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_japanese_ricefish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_japanese_ricefish/compare_total_rep5.png
cd ../../../
echo "Zebrafish"
cd 1e0o_zebrafish/charmm-gui-8224620512/
cd namdrep1/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_zebrafish/compare_total_rep1.png
cd ../namdrep2/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_zebrafish/compare_total_rep2.png
cd ../namdrep3/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_zebrafish/compare_total_rep3.png
cd ../namdrep4/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_zebrafish/compare_total_rep4.png
cd ../namdrep5/
cp /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/analyze_mmpbsa.py .
python3 analyze_mmpbsa.py summary iface1_run/FINAL_RESULTS_iface1.dat
python3 analyze_mmpbsa.py compare \
    --label iface1 iface1_run/FINAL_RESULTS_iface1.dat \
    --label iface2 iface2_run/FINAL_RESULTS_iface2.dat \
    --label whole  whole_run/FINAL_RESULTS_whole.dat \
    -o /Users/alexanderwang/Desktop/md_sims/all_metrics_analysis/mmpbsa_numbers/1e0o_zebrafish/compare_total_rep5.png
cd ../../../
