import pandas as pd
import re
from collections import defaultdict

# 读取defense_islands.txt文件
with open('5455MAG_defense_islands.txt', 'r') as file:
    lines = file.readlines()

# 初始化变量
defense_islands = []
current_island = None
all_defense_types = set()

# 解析文件内容
for line in lines:
    line = line.strip()
    if line.startswith("Defense Island"):
        if current_island:
            defense_islands.append(current_island)
        current_island = {
            'name': line.split(':')[0],
            'contig': '',
            'gene_count': 0,
            'defense_types': defaultdict(int),
            'first_record_num': None,
            'last_record_num': None,
            'gene_set': set() #初始化基因集合
        }
    elif line:
        parts = line.split('\t')
        contig = parts[0]
        defense_type = parts[1]
        # 从parts[5]读取以逗号分隔的字符串
        current_genes = set(parts[6].split(',')) if len(parts) > 6 else set()  # 添加检查以防止索引错误
        # 将当前行的基因添加到当前岛屿的基因集中
        current_island['gene_set'].update(current_genes)  # 更新基因集合
        
        first_record_num = int(parts[4].split('_')[-1])
        last_record_num = int(parts[5].split('_')[-1])
        
        current_island['contig'] = contig
        current_island['gene_count'] = len(current_island['gene_set'])  # 计算非冗余基因数量
        current_island['defense_types'][defense_type] += 1  # 统计当前防御类型的数量
        all_defense_types.add(defense_type)
        
        if current_island['first_record_num'] is None:
            current_island['first_record_num'] = first_record_num
        # 仅在新的last_record_num更大时更新
        if current_island['last_record_num'] is None or last_record_num > current_island['last_record_num']:
            current_island['last_record_num'] = last_record_num

# 添加最后一个岛屿
if current_island:
    defense_islands.append(current_island)

# 创建DataFrame
data = []
for island in defense_islands:
    row = [
        island['name'],
        island['contig'],
        island['gene_count'],
        ','.join(island['defense_types'].keys()),
        island['first_record_num'],
        island['last_record_num']
    ]
    for defense_type in all_defense_types:
        row.append(island['defense_types'][defense_type])
    data.append(row)

columns = [
    'Defense Island Name',
    'Contig',
    'Total Gene Count',
    'Defense Types',
    'First Record Num',
    'Last Record Num'
] + list(all_defense_types)

df = pd.DataFrame(data, columns=columns)

# 保存为CSV文件
df.to_csv('5455MAG_defense_islands_summary.csv', sep='\t', index=False)

# 打印DataFrame
print(df)