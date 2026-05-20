import pandas as pd

# 读取CSV文件
file_path = '5455MAG_defense_islands_summary.csv'
df = pd.read_csv(file_path, sep='\t')

# 处理Contig列
def process_contig(row):
    # 删除'bin.'之后的部分，保留'bin.'后紧跟的数字和下划线
    contig_parts = row['Contig'].split('bin.')
    if len(contig_parts) > 1:
        contig = contig_parts[0] + 'bin.' + contig_parts[1].split('_')[0] + '_'
    else:
        contig = row['Contig']  # 如果没有'bin.'，则保持原样
    
    # 生成从First Record Num到Last Record Num的数字字符串
    record_nums = list(range(row['First Record Num'], row['Last Record Num'] + 1))
    new_contigs = [f"{contig}{num}" for num in record_nums]  # 每个数字单独添加到contig后
    return new_contigs  # 返回生成的contig列表

# 应用处理函数并展开结果
all_contigs = []
for _, row in df.iterrows():
    all_contigs.extend(process_contig(row))  # 将每行生成的contig添加到总列表中

# 输出到txt文件，每个contig一行
output_file_path = '01all_DIorf.txt'
with open(output_file_path, 'w') as f:
    for contig in all_contigs:
        f.write(f"{contig}\n")

print("处理完成，结果已输出到 01all_DIorf.txt。")
