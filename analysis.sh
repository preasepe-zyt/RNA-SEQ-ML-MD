source /usr/local/gromacs/bin/GMXRC
# 获取用户输入
read -p "请输入一个值: " user_input

echo 4 4 | gmx rms -s  ../md_0_1.tpr -f  ../md_0_1.xtc -o rmsd_$user_input.xvg -fit rot+trans -xvg none
echo 1 | gmx rmsf -s ../md_0_1.tpr -f ../md_0_1.xtc  -o rmsf_$user_input.xvg -fit -xvg none
echo 1 | gmx gyrate -s ../md_0_1.tpr -f ../md_0_1.xtc -o gyrate_$user_input.xvg -xvg none
echo 1 | gmx sasa -s ../md_0_1.tpr -f ../md_0_1.xtc -o area_$user_input.xvg -or resarea.xvg -oa atomarea.xvg -xvg none

conda activate gmxMMPBSA
gmx mdrun -rerun ../md_0_1.trr -s ../md_0_1.tpr -o rerun.trr -c rerun.gro
#提取xtc文件
gmx  trjconv -f rerun.trr -o  mdreturn.xtc
#矫正
echo 1 0 | gmx trjconv -s ../md_0_1.tpr -f mdreturn.xtc -ur rect -pbc mol -center -o md_center.xtc
#蛋白质二级结构
echo 1 | gmx do_dssp -f  md_center.xtc -s ../md_0_1.tpr -n ../index.ndx -tu ns -o ss.xpm -sc scount_$user_input.xvg  -xvg none
# 设置变量
START_FRAME=1
END_FRAME=5000
INTERVAL=5
VERBOSE=2
FORCEFIELDS="oldff/leaprc.ff99SB,leaprc.gaff"
IGB=5
SALTCON=0.150

# 创建 mmpbsa.in 文件
cat <<EOF > mmpbsa.in
&general
startframe=1, endframe=5000, interval=5, verbose=2, 
forcefields="oldff/leaprc.ff99SB,leaprc.gaff"
/
&gb
igb=5, saltcon=0.150
/
&decomp
idecomp=2, dec_verbose=3,
print_res="within 4"
/
EOF

# 运行 gmx_MMPBSA
mpirun -np 10 gmx_MMPBSA MPI -O -i mmpbsa.in -cs ../md_0_1.tpr -ci  ../index.ndx -cg 1 13 -ct md_center.xtc -cp ../topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv
conda deactivate




