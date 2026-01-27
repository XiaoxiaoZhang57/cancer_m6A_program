import numpy as np
from collections import Counter

# 定义标准氨基酸列表 (20种标准氨基酸)
STANDARD_AAS = list("ACDEFGHIKLMNPQRSTVWY")

def read_alignment(file_path):
    """读取CLUSTAL格式的比对文件"""
    sequences = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # 跳过前两行（标题和空行）
    for line in lines[2:]:
        if line.strip():  # 跳过空行
            parts = line.split()
            if len(parts) > 1:
                # 取序列部分（忽略序列名）
                seq = parts[1].upper()
                sequences.append(seq)
    
    return sequences

def calculate_conservation_score(column, p_e=0.05):
    """
    计算单个位点的保守性分数
    M(i) = Σ [p_i(x) - p_e]² / p_e
    """
    # 过滤掉gap和非标准氨基酸
    valid_aas = [aa for aa in column if aa in STANDARD_AAS]
    
    if not valid_aas:
        return 0.0  # 如果没有有效氨基酸，返回0
    
    n = len(valid_aas)
    counter = Counter(valid_aas)
    
    score = 0.0
    for aa in STANDARD_AAS:
        p_i = counter.get(aa, 0) / n  # 观测概率
        deviation = p_i - p_e
        score += (deviation * deviation) / p_e
    
    return score

def main():
    # 读取比对文件
    sequences = read_alignment("ENST00000375486.9.long24hg")
    
    # 获取比对长度（以第一条序列为准）
    alignment_length = len(sequences[0])
    
    # 计算每个位点的保守性分数
    conservation_scores = []
    for i in range(alignment_length):
        # 获取当前位点的所有氨基酸
        column = [seq[i] for seq in sequences if i < len(seq)]
        score = calculate_conservation_score(column)
        conservation_scores.append(score)
    
    # 输出结果
    print("Position\tConservation Score")
    for i, score in enumerate(conservation_scores, 1):
        print(f"{i}\t\t{score:.4f}")
    
    # 可选：保存结果到文件
    with open("conservation_scores.txt", "w") as f:
        f.write("Position\tConservation Score\n")
        for i, score in enumerate(conservation_scores, 1):
            f.write(f"{i}\t\t{score:.4f}\n")

if __name__ == "__main__":
    main()
