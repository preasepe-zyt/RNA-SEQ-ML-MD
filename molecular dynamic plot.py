import matplotlib.pyplot as plt
import matplotlib
import numpy
import sys
sys.path.append(r"C:\Users\79403\Desktop\fsdownload\zhao")
from mean import mean
import matplotlib.pyplot as plt
import os
from str_split import str_split
import pandas as pd

#matplotlib.use("TKAgg")
#字体
#matplotlib.rcParams['font.family'] = 'Arial'
plt.rcParams["font.family"] = "serif"
#plt.rcParams["font.serif"] = "Times New Roman"
path = r"C:\Users\79403\Desktop\fsdownload\zhao"
# 定义颜色列表
colors = ['blue', 'green', 'red', 'purple', 'orange', "olive"]
#RMSD
list_rmsd = os.listdir(path)
rmsd_name = [i for i in list_rmsd if "rmsd" in i]
rmsd_legend = [str_split(i,"_",".") for i in rmsd_name]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)
for index, i in enumerate(rmsd_name):
    t, rmsd = numpy.loadtxt(path + "\\" + i, unpack=True)
    t_m, rmsd_m = mean(rmsd)
    ax.plot(t/1000, rmsd, linestyle="-", alpha=0.2, color=colors[index])
    ax.plot(t_m/100, rmsd_m, linestyle="-", color=colors[index], label=rmsd_legend[index]) #
ax.set_xlabel("$t/ns$", fontsize=20)
ax.set_ylabel(r"RMSD (nm)", fontsize=20) #r"C$_\alpha$ RMSD (nm)"
ax.legend(loc='upper right', frameon=False, fontsize=13)
# 设置坐标轴刻度的大小
ax.tick_params(axis='both', labelsize=20)
# 设置坐标轴边框粗细
for spine in ax.spines.values():
    spine.set_linewidth(2)
ax.set_yticks(numpy.arange(0, 1.2, 0.1))
plt.show()
fig.savefig(path+"\\"+"rmsd.png", dpi=300)

#rmsf
import matplotlib.pyplot as plt
import numpy as np
list_rmsd = os.listdir(path)
rmsf_name = [i for i in list_rmsd if "rmsf" in i]
rmsf_legend = [str_split(i,"_",".") for i in rmsf_name]
residue = ("LEU289", "VAL244", "ALA454", "LEU67",
           "VAL27", "VAL111", "HIS131", "ILE127",
           "ARG102", "GLN142")
residue_lo = [list(range(3699, 3703)), list(range(2035, 2041)), list(range(2401, 2416)),  list(range(469, 476)),
              list(range(196, 202)), list(range(900, 906)), list(range(1064, 1079)), list(range(1030, 1073))
              , list(range(825, 835)), list(range(1167, 1175))]
style = {'color': 'black', 'fontsize': 11, 'fontweight': 'bold'}
resisue_h = [0.4, 0.42, 0.55, 0.4, 0.53, 0.68, 0.52, 0.3, 0.6, 0.73]

fig, axs = plt.subplots(2, 3, figsize=(12, 8))
fig.subplots_adjust(wspace=0.4, hspace=0.4)
for index, i in enumerate(rmsf_name):
    # 计算子图的索引位置
    row = index // 3 #整数除法
    col = index % 3 #求余数

    # 加载数据，这里使用随机数据代替
    resid, rmsf = np.loadtxt(path+"\\"+i, unpack=True)

    # 在相应的子图中绘制数据
    axs[row, col].plot(resid, rmsf, linestyle="-", linewidth=0.5, color=colors[index])
    axs[row, col].set_xlabel("Residue number", fontsize=20)
    axs[row, col].set_title(rmsf_legend[index])
    axs[row, col].set_ylabel("RMSF (nm)", fontsize=20)
    axs[row, col].tick_params(axis='both', labelsize=10)
    axs[row, col].spines['top'].set_linewidth(2)
    axs[row, col].spines['right'].set_linewidth(2)
    axs[row, col].legend(loc='upper right', frameon=False)
plt.show()
fig.savefig(path+"\\"+"rmsf.png", dpi=300)


#distance
import matplotlib.pyplot as plt
import numpy

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

for index, i in enumerate(rmsd_name):
   t, d, a, b = numpy.loadtxt(path+"\\"+i, unpack=True)
   ax.plot(t / 1000, d, linestyle="-")

ax.set_xlabel("$t/ns$", fontsize=20)
ax.set_ylabel("Distance (nm)", fontsize=20)
ax.tick_params(axis='both', labelsize=20)
for spine in ax.spines.values():
    spine.set_linewidth(2)
plt.show()
fig.savefig("dis.png", dpi=300)

#Radius of gyration
import matplotlib.pyplot as plt
import numpy
list = os.listdir(path)
gyrate_name = [i for i in list if "gyrate" in i]
gyrate_legend = [str_split(i,"_",".") for i in gyrate_name]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

for index, i in enumerate(gyrate_name):
    t, data, x, y, z = numpy.loadtxt(path+"\\"+i, unpack=True)
    t_m, data_m = mean(data)
    ax.plot(t / 1000, data, linestyle="-", color=colors[index], alpha=0.2)
    ax.plot(t_m / 100, data_m, linestyle="-", color=colors[index], label=gyrate_legend[index])


ax.set_xlabel("$t/ns$", fontsize=20)#$s$ 可以斜体
ax.set_ylabel(r"RoG (nm)", fontsize=20)# _\mathrm{hb} 可以把字体下标
ax.legend(loc='upper right', frameon=False, fontsize=13)
ax.set_yticks(numpy.arange(2.40, 3.80, 0.10))
# 设置坐标轴刻度的大小
ax.tick_params(axis='both', labelsize=20)
# 设置坐标轴边框粗细
for spine in ax.spines.values():
    spine.set_linewidth(2)
plt.show()
fig.savefig(path+"\\"+"rgyr.png", dpi=300)

#fig.savefig("rgyr.svg")
#fig.savefig("rgyr.pdf")


#hbond
import matplotlib.pyplot as plt
import numpy
list = os.listdir(path)
hbnum_name = [i for i in list if "hbnum" in i]
hbnum_legend = [str_split(i,"_",".") for i in hbnum_name]
fig = plt.figure(figsize=(8, 4))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)
for index, i in enumerate(hbnum_name):
    t, x, hbond = numpy.loadtxt(path + "\\" + i, unpack=True)
    ax.bar(t/1000, x, alpha=0.5, label=hbnum_legend[index])
#plt.legend(loc='upper left', frameon=False)
ax.set_xlabel("$t/ns$", fontsize=20)#$s$ 可以斜体
ax.set_ylabel("Number of Hbond", fontsize=20)# _\mathrm{hb} 可以把字体下标
ax.set_yticks(numpy.arange(0, 8, 1))
# 设置坐标轴刻度的大小
ax.tick_params(axis='both', labelsize=20)
ax.legend(loc='upper right', frameon=False, fontsize=13)
# 设置坐标轴边框粗细
for spine in ax.spines.values():
    spine.set_linewidth(2)
plt.show()
fig.savefig(path+"\\"+"hbond.png", dpi=300)
#
'tab:blue'
'tab:orange'
'tab:green'
'tab:red'
'tab:purple'
'tab:brown'
'tab:pink'
'tab:gray'
'tab:olive'
'tab:cyan'

#绘制二级结构
import numpy as np
import matplotlib.pyplot as plt
# color_map = ["总结构 (Structure)", "区卷 (Coil)", "β-折叠 (B-Sheet)","β-桥 (B-Bridge)", "弯曲 (Bend)",
#               "转弯 (Turn)", "α-螺旋 (A-Helix)", "5-螺旋 (5-Helix)", "分链器 (3-Helix)"]
color_map = ["Structure", "Coil", "B-Sheet", "B-Bridge", "Bend",
              "Turn", "A-Helix", "5-Helix", "3-Helix"]
# 替换为你的 DSSP 输出文件路径
fig, axs = plt.subplots(2, 3, figsize=(12, 12))
fig.subplots_adjust(wspace=0.45, hspace=0.6)
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定中文字体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

for index, i in enumerate(group):
    # 计算子图的索引位置
    row = index // 3  # 整数除法
    col = index % 3  # 求余数
    # 加载数据，这里使用随机数据代替
    data = np.loadtxt(path + "\\" + "scount" + "_" + i + ".xvg", unpack=True)
    for a in range(1, 9):
        axs[row, col].plot(data[0], data[a], linestyle="-", label=color_map[a-1], linewidth=0.5)
    axs[row, col].set_xlabel("t/ns", fontsize=20)
    axs[row, col].set_title(i)
    axs[row, col].set_ylabel("Residue", fontsize=20)
    axs[row, col].tick_params(axis='both', labelsize=9)
    axs[row, col].spines['top'].set_linewidth(2)
    axs[row, col].spines['right'].set_linewidth(2)
    axs[row, col].legend(loc='upper center', bbox_to_anchor=(0.5, 1.4), frameon=False, title="Second Structure", ncol=3)
plt.show()
fig.savefig(path+"\\"+"second_structure.png", dpi=300)


#溶剂可及面积 SASA
list = os.listdir(path)
area_name = [i for i in list if "area" in i]
area_legend = [str_split(i,"_",".") for i in gyrate_name]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

for index, i in enumerate(area_name):
    t, sasa = numpy.loadtxt(path+"\\"+i, unpack=True)
    t_m, sasa_m = mean(sasa)
    ax.plot(t/1000, sasa, linestyle="-", color=colors[index], alpha=0.2)
    ax.plot(t_m / 100, sasa_m, linestyle="-", color=colors[index], label=area_legend[index])

ax.set_xlabel("$t/ns$", fontsize=20)
ax.set_ylabel(r"SASA (nm²)", fontsize=20) #r"C$_\alpha$ RMSD (nm)"
ax.legend(loc='upper right',frameon=False, fontsize=13)
# 设置坐标轴刻度的大小
ax.tick_params(axis='both', labelsize=20)
# 设置坐标轴边框粗细
for spine in ax.spines.values():
    spine.set_linewidth(2)
ax.set_yticks(numpy.arange(260, 400, 10))
plt.show()
fig.savefig(path+"\\"+"area.png", dpi=300)

#结合自由能
list = os.listdir(path)
MMPBSA_name = [i for i in list if "MMPBSA" in i]
#MMPBSA_legend = [str_split(i,"_",".") for i in hbnum_name]
data = pd.read_csv(path+"\\"+MMPBSA_name[0],sep=",")
col = ['Frame #','VDWAALS', 'EEL','EGB', 'ESURF', 'TOTAL']
data = data[col]
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

for index, i in enumerate(data.columns.drop('Frame #')):
    ax.plot((data.iloc[:,0]*50)/100, data[i], linestyle="-", color=colors[index], label=i)

ax.set_xlabel("$t/ns$", fontsize=20)
ax.set_ylabel(r"Binding Energy (kcal/mol)", fontsize=20) #r"C$_\alpha$ RMSD (nm)"
ax.legend(loc='upper center', frameon=False, fontsize=13, ncol=5, bbox_to_anchor=(0.5, 1.1))
# 设置坐标轴刻度的大小
ax.tick_params(axis='both', labelsize=20)
# 设置坐标轴边框粗细
for spine in ax.spines.values():
    spine.set_linewidth(2)
#ax.set_yticks(numpy.arange(-60, 60, 10))
plt.show()
fig.savefig(path+"\\"+"MMPBSA.png", dpi=300)