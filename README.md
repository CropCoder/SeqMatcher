# SeqMatcher

高性能多线程 DNA 序列引物匹配与文库变体计数工具。

将海量测序序列与已知引物库进行匹配，统计各引物的序列覆盖度，并量化文库变体在匹配序列中的出现频次。全链路 Rust 实现，Rayon 多线程并行，核心热路径用编译期查找表优化。

## 安装

```bash
git clone https://github.com/CropCoder/SeqMatcher.git
cd SeqMatcher
cargo build --release
```

二进制文件位于 `target/release/seq_matcher`。

## 快速开始

```bash
./target/release/seq_matcher \
  --primer-csv primers.csv \
  --library-csv library.csv \
  --seq a_11:data/11_seq.txt \
  --output-dir output
```

### 输入文件格式

**引物 CSV**（`--primer-csv`）— 三列：引物 ID、正向序列、反向序列。

| primer_id | forward_seq | reverse_seq |
|-----------|-------------|-------------|
| P001      | ATCGGTACC   | GCTATAGCA   |
| P002      | TGCACTGAC   | CGTACGATG   |

**文库 CSV**（`--library-csv`）— 含一列变体序列，列名可配置。

| variant_id | single_degenerate_library_expanded_reference |
|------------|---------------------------------------------|
| V001       | ATCGNNNTCGA                                 |
| V002       | GCTANNNGGCTA                                 |

**序列文件**（`--seq`）— 每行一条序列的纯文本文件。

### 输出文件

每次运行对每个 `--seq` 输入生成两个 CSV：

- `{LABEL}_seq_matched_primers_count.csv` — 各引物匹配到的序列总数
- `{LABEL}_seq_matched_library_variant_count.csv` — 原始文库表 + 每个引物的变体命中计数列

## CLI 参考

```
Usage: seq_matcher [OPTIONS] --primer-csv <PRIMER_CSV> --library-csv <LIBRARY_CSV>

Options:
  -p, --primer-csv <PRIMER_CSV>              引物 CSV 文件路径 (列: id, forward_seq, reverse_seq)
  -l, --library-csv <LIBRARY_CSV>            文库 CSV 文件路径
      --library-seq-col <LIBRARY_SEQ_COL>    文库 CSV 中序列所在列名
                                              [default: single_degenerate_library_expanded_reference]
  -s, --seq <SEQ_FILES>                      序列文件: 格式为 LABEL:PATH (如 a_11:data/11_seq.txt)
                                              可多次指定以批量处理
  -o, --output-dir <OUTPUT_DIR>              输出目录 [default: output]
  -c, --chunk-size <CHUNK_SIZE>              并行处理块大小 (条/批) [default: 10000]
  -t, --threads <THREADS>                    线程数 (默认使用全部 CPU 核心)
  -h, --help                                 打印帮助信息
  -V, --version                              打印版本号
```

## 算法

1. 加载引物表与文库表，预计算所有序列的反向互补。
2. 目标序列按可配块大小分片，Rayon 多线程并行处理。
3. 每条序列**首次命中**引物即停止检索（first-match-wins）。
4. 匹配到引物后，扫描全部文库变体，检查其原始序列或反向互补是否**包含于**目标序列中。
5. 各线程内部无锁统计，块处理后合并至全局结果。
6. 批量输出两个 CSV 文件。

### 性能

- 反向互补使用编译期 `const` 128 字节查找表，单周期映射。
- 共享数据（引物、文库）用 `Arc` 零拷贝跨线程传递。
- 输出使用 `BufWriter` 批量写入，避免逐行 IO。

## 依赖

| crate  | 用途           |
|--------|--------------|
| clap   | CLI 参数解析    |
| csv    | CSV 读写      |
| rayon  | 数据并行        |
| serde  | 序列化 (derive) |
| anyhow | 错误处理        |

## License

MIT
