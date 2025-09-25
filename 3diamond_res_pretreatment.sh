
#根据ORF提取MAG
cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/1ARG
awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_ARG_out > diamond_ARG_out2
awk 'BEGIN{OFS="\t"} {print $0, "ARG"}' diamond_ARG_out2 > diamond_ARG_out3
python3 /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py diamond_ARG_out3

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/2VFG
awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_VFG_out > diamond_VFG_out2
awk 'BEGIN{OFS="\t"} {print $0, "VFG"}' diamond_VFG_out2 > diamond_VFG_out3
python3 /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py diamond_VFG_out3

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/3MRG
awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_MRG_out > diamond_MRG_out2
awk 'BEGIN{OFS="\t"} {print $0, "MRG"}' diamond_MRG_out2 > diamond_MRG_out3
python3 /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py diamond_MRG_out3

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4mibPOP
awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    # 在原行末尾追加提取的内容
    print $0, extracted
}' diamond_POPDeg_out > diamond_POPDeg_out2
awk 'BEGIN{OFS="\t"} {print $0, "POPDeg"}' diamond_POPDeg_out2 > diamond_POPDeg_out3
python3 /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py diamond_POPDeg_out3

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/0ARG_RGI

cut -f1,15,16,17,28 rgiargs.txt > rgiargs_2.txt

awk 'BEGIN{OFS="\t"} {
    # 从第一列提取 "bin.数字" 之前的所有内容
    extracted = "NA"  # 默认值
    if (match($1, /.*bin\.[0-9]+/)) {
        extracted = substr($1, 1, RLENGTH)
    }
    $1 = extracted
	print $0
}' rgiargs_2.txt > rgiargs_3.txt

python3 /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/4turn_to_table.py rgiargs_2.txt

