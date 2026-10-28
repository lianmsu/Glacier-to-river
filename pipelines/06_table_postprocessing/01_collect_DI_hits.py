import pandas as pd
import glob
import re

# 定义提取contig_id的函数
def extract_contig_id(sys_id):
    match = re.match(r'(.+-bin.\d+\-\d+)', sys_id)
    if match:
        return match.group(1)
    return sys_id

# 初始化变量
defense_islands = []

# 获取所有TSV文件的列表
tsv_files = glob.glob('/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/DIs/*.tsv')

# 遍历每个TSV文件
for tsv_file in tsv_files:
#    print(f"处理文件: {tsv_file}")
    
    # 读取TSV文件
    df = pd.read_csv(tsv_file, sep='\t')
    
    #只保留activity为“Defense”的行
    #df = df[df['activity'] == 'Defense'].copy()  # 使用 .copy() 避免 SettingWithCopyWarning

    # 提取contig ID
    df['contig_id'] = df['sys_id'].apply(extract_contig_id)
    
    # 打印前几行数据以检查
#    print("DataFrame 前几行数据:")
#    print(df.head())
    
    # 按contig_id分组
    grouped = df.groupby('contig_id')
    
    for contig_id, group in grouped:
#        print(f"\n处理 contig_id: {contig_id}")
        # 提取sys_beg末尾的数字并排序
        group = group.assign(sys_beg_num=group['sys_beg'].apply(lambda x: int(x.split('_')[-1])))
        group = group.sort_values(by='sys_beg_num')
        current_island = []
        previous_gene_end = None

        # 遍历每一行
        for index, row in group.iterrows():
            current_gene_start = int(row['sys_beg'].split('_')[-1])
            current_gene_end = int(row['sys_end'].split('_')[-1])
            
            if previous_gene_end is None or current_gene_start - previous_gene_end <= 10:
                current_island.append(row)
            else:
                if sum(gene['genes_count'] for gene in current_island) >= 5 and len(set([gene['subtype'] for gene in current_island])) >= 3:
                    defense_islands.append(current_island)
#                    print(f"找到防御岛屿: {current_island}")
                current_island = [row]
            
            previous_gene_end = current_gene_end

        # 检查最后一个岛屿
        if sum(gene['genes_count'] for gene in current_island) >= 5 and len(set([gene['subtype'] for gene in current_island])) >= 3:
            defense_islands.append(current_island)
#            print(f"找到防御岛屿: {current_island}")

# 输出结果到txt文件
with open('5455MAG_defense_islands.txt', 'w') as f:
    for i, island in enumerate(defense_islands):
        f.write(f"Defense Island {i+1}:\n")
        for gene in island:
            f.write('\t'.join(map(str, gene.values)) + '\n')
        f.write('\n')

# 打印防御岛屿数量
print(f"\n找到的防御岛屿数量: {len(defense_islands)}")