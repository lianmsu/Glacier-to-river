#################1###############################
import pandas as pd

# 读取CSV文件
file_path = 'defense_islands_summary4.csv'
df = pd.read_csv(file_path, sep='\t')

# 处理Contig列
def process_contig(row):
    # 删除'cov_'之后的部分，保留'cov_'后紧跟的数字和下划线
    contig_parts = row['Contig'].split('cov_')
    if len(contig_parts) > 1:
        contig = contig_parts[0] + 'cov_' + contig_parts[1].split('_')[0] + '_'
    else:
        contig = row['Contig']  # 如果没有'cov_'，则保持原样
    
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

###########################2##########################################

#1. 生成所有DIs所在的contig的所有orf id
python3 /home/lipengwei2024phd/03GW/02ANI99MAGs_defense/DIs_unknowngene/01all_DI_orf.py
#2. 提取所有DI对应的DGs的id
cat /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/DI_othergenes/5455MAG_defense_islands.txt | cut -f7 | grep -v ":" | grep -v '^ *$' | tr ',' '\n' | sort | uniq > 02DI_contig_DGs.txt
#注意：defensefinder在识别防御系统时会把一个DGs重复算到不同的DSs

##################################3##################################

#1. 生成所有DIs所在的contig的所有orf id
python3 /home/lipengwei2024phd/03GW/02ANI99MAGs_defense/DIs_unknowngene/01all_DI_orf.py

#2. 提取所有DI对应的DGs的id
cat /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/DI_othergenes/5455MAG_defense_islands.txt | cut -f7 | grep -v ":" | grep -v '^ *$' | tr ',' '\n' | sort | uniq > 02DI_contig_DGs.txt
#注意：defensefinder在识别防御系统时会把一个DGs重复算到不同的DSs

#3. 从01all_DIorf.txt中删除与02DI_contig_DGs.txt共有的行
# 首先，提取并排序两个文件的每一行
sort 01all_DIorf.txt | uniq > 03all_DIorf_sorted.txt
sort 02DI_contig_DGs.txt | uniq > 03DI_contig_DGs_sorted.txt

comm -23 03all_DIorf_sorted.txt 02DI_contig_DGs.txt > 03filtered_01all_DIorf.txt

# 找出两个排序文件共有的行
#grep -Ff 03DI_contig_DGs_sorted.txt 03all_DIorf_sorted.txt > 03common_lines.txt
# 从原始文件01all_DIorf.txt中删除这些共有行
#grep -v -Ff 03common_lines.txt 01all_DIorf.txt > 03all_DIorf_filtered.txt

#4. 根据03all_DIorf_filtered.txt提取orf序列
seqkit grep -f 03filtered_01all_DIorf.txt /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG.prt_defensefinder.prt > 04non_DG_orf.faa

#根据文件03filtered_01all_DIorf.txt提取对应的gff记录
#!/bin/bash
cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/DI_othergenes
# 输入文件路径
gene_list="/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/DI_othergenes/03filtered_01all_DIorf.txt"
gff_file="/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG.prt_defensefinder.gff"
output_file="filtered_DI_genes.gff"

# 步骤1: 创建基因ID列表的临时文件（去除可能的空行）
grep -v "^$" "$gene_list" | sort -u > gene_ids.tmp

# 步骤2: 从GFF文件中提取匹配的记录
echo "Extracting matching records from GFF file..."
awk -F"\t" 'BEGIN{OFS="\t"} 
NR==FNR {
    # 读取基因ID列表到数组
    genes[$1]=1;
    next;
}
/^#/ {
    # 保留GFF头信息
    print;
    next;
}
{
    # 提取ID属性
    split($9, attrs, ";");
    gene_id="";
    for(i in attrs) {
        if(attrs[i] ~ /^ID=/) {
            split(attrs[i], id, "=");
            gene_id=id[2];
            break;
        }
    }
    if(gene_id in genes) {
        print;
    }
}' gene_ids.tmp "$gff_file" > "$output_file"

# 清理临时文件
rm gene_ids.tmp

echo "Extraction complete! Results saved to: $output_file"
echo "Number of matched records: $(grep -v "^#" "$output_file" | wc -l)"


#5. 在线eggnog-mapper注释

##################################

