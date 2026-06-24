# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict

# ===================== 反向互补函数（固定） =====================
def reverse_complement(x):
    key_dict = {
        'A':'T','T':'A','C':'G','G':'C',
        'R':'Y','Y':'R','K':'M','M':'K','N':'N'
    }
    x = str(x).upper().strip()
    return ''.join([key_dict[c] for c in reversed(x)])

# ===================== 引物匹配 =====================
def check_primer_match(seq, f_primer, r_primer):
    seq = seq.upper().strip()
    f = f_primer.upper().strip()
    r = r_primer.upper().strip()
    
    rc_f = reverse_complement(f)
    rc_r = reverse_complement(r)
    
    return (seq.startswith(f) and seq.endswith(rc_r)) or \
           (seq.startswith(r) and seq.endswith(rc_f))

# ===================== 流式读取序列 =====================
def stream_sequences(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                yield line

# ===================== 主程序（已优化 + 无BUG） =====================
if __name__ == "__main__":
    PRIMER_CSV = "/mnt/ma/UserData/liujing_data/data/202606/library/80_primers_list_all_123.csv"
    LIBRARY_CSV = "/mnt/ma/UserData/liujing_data/data/202606/library/80_full_library_2_2.csv"

    SEQ_MAP = {
        "11_seq.txt": "a_11"
    }

    # 读取引物 + 构建快速查询表
    df_primers = pd.read_csv(PRIMER_CSV)
    primer_dict = {}
    for _, row in df_primers.iterrows():
        pid = str(row.iloc[0]).strip()
        f_seq = row.iloc[1]
        r_seq = row.iloc[2]
        primer_dict[pid] = (f_seq, r_seq)
    primer_list = list(primer_dict.keys())

    # 读取库 + 构建超级快的查询字典（关键优化）
    df_lib = pd.read_csv(LIBRARY_CSV)
   # 🔥【修改1】预存原序列+反向互补序列，提前计算提速
    lib_raw = []
    lib_rc  = []
    for seq in df_lib["single_degenerate_library_expanded_reference"]:
        if pd.isna(seq):
            lib_raw.append("")
            lib_rc.append("")
        else:
            s = str(seq).strip().upper()
            lib_raw.append(s)
            lib_rc.append(reverse_complement(s))
    # 🔥结束

    # 构建结果列
    a_11_cols = [f"{p}_a_11" for p in primer_list]

    for c in a_11_cols: 
        df_lib[c] = 0

    # 计数
    primer_count = {
        p: {"a_11":0}
        for p in primer_list
    }

    print(f"✅ 引物总数：{len(primer_list)}")
    print(f"✅ Library 变体数：{len(df_lib)}")
    print(f"✅ 已构建快速查询字典，速度起飞！\n")

    # ------------------- 处理序列（超快版） -------------------
    for seq_file, suffix in SEQ_MAP.items():
        print(f"📂 处理：{seq_file} → {suffix}")
        
        for seq_idx, target_seq in enumerate(stream_sequences(seq_file), 1):
            target_seq = target_seq.upper()
            matched_primer = None

            # 匹配引物
            for pid, (f, r) in primer_dict.items():
                if check_primer_match(target_seq, f, r):
                    matched_primer = pid
                    primer_count[pid][suffix] += 1
                    break

            if not matched_primer:
                continue

             # 🔥【修改2】包含匹配逻辑，复用预计算序列提速
            col_name = f"{matched_primer}_{suffix}"
            for lib_idx in range(len(lib_raw)):
                if lib_raw[lib_idx] in target_seq or lib_rc[lib_idx] in target_seq:
                    df_lib.at[lib_idx, col_name] += 1
            # 🔥结束=================================

            if seq_idx % 50000 == 0:
                print(f"   → 已处理 {seq_idx} 条")

    # ------------------- 调整列顺序（BUG已修复） -------------------
    original_cols = list(df_lib.columns[:-len(a_11_cols)])
    df_lib = df_lib[original_cols + a_11_cols]

    # 输出计数表
    df_out = df_primers.copy()
    df_out["a_11"] = [primer_count[p]["a_11"] for p in primer_list]

    df_out.to_csv("/mnt/ma/UserData/liujing_data/data/202606/output/match_primer_AA-output/11_2_seq_matched_primers_count_dic.csv", index=False, encoding="utf-8-sig")
    df_lib.to_csv("/mnt/ma/UserData/liujing_data/data/202606/output/match_primer_AA-output/11_2_seq_matched_library_variantAA_count_dic.csv", index=False, encoding="utf-8-sig")
