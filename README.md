# SeqMatcher

高性能多线程 DNA 序列引物匹配与文库变体计数工具。

将海量测序序列与已知引物库进行匹配，统计各引物的序列覆盖度，并量化文库变体在匹配序列中的出现频次。Rust 全链路实现，Rayon 多线程并行 + Aho-Corasick 多模式匹配，核心热路径用编译期查找表优化。

## 安装

```bash
git clone https://github.com/CropCoder/SeqMatcher.git
cd SeqMatcher
cargo build --release
```

二进制文件位于 `target/release/seq_matcher`。

### macOS → Linux 交叉编译

```bash
# 前提: brew install musl-cross && rustup target add x86_64-unknown-linux-musl
./build-linux.sh
```

生成 `target/x86_64-unknown-linux-musl/release/seq_matcher`，ELF 64-bit 静态链接，可直接拷贝到 Linux 运行。

## 快速开始

```bash
./target/release/seq_matcher \
  --primer-csv primers.csv \
  --library-csv library.csv \
  --seq a_11:data/11_seq.txt \
  --output-dir output
```

### 运行示例

```
  Loading primers from: primers_list_all.csv
  Loaded 32 primers
  Loading library from: 80_full_library_2_12.csv
  Loaded 8192 library variants
  Processing: data/11_seq.txt -> a_11
    总条数: 1500000  |  chunk: 10000  |  引物: 32  |  变体: 8192  |  AC patterns: 16384
    [████████████████████░░░░░░░░░░░░░░░░]  55.0%  825000/1500000  45230 seq/s  ETA: 15s
    完成: 1500000 条序列, 耗时 33.2s, 速度 45181 seq/s
  All done.
```

### 输入文件格式

**引物 CSV**（`--primer-csv`）— 任意列数，工具保留所有原始列。默认取前三列为引物 ID、正向序列、反向序列。

| primer_id | forward_seq | reverse_seq | 其他列...          |
|-----------|-------------|-------------|-------------------|
| P001      | ATCGGTACC   | GCTATAGCA   | （保留，原样输出）    |
| P002      | TGCACTGAC   | CGTACGATG   | （保留，原样输出）    |

**文库 CSV**（`--library-csv`）— 含一列变体序列，列名通过 `--library-seq-col` 配置。

| variant_id | single_degenerate_library_expanded_reference | 其他列...       |
|------------|---------------------------------------------|-----------------|
| V001       | ATCGNNNTCGA                                 | （保留，原样输出） |

**序列文件**（`--seq`）— 每行一条 DNA 序列的纯文本文件。

### 输出文件

每次运行对每个 `--seq` 输入生成两个 CSV：

- `{LABEL}_seq_matched_primers_count.csv` — 原始引物表所有列 + 一列 `count_{LABEL}`（每条引物匹配到的序列总数）
- `{LABEL}_seq_matched_library_variant_count.csv` — 原始文库表所有列 + 每个引物的变体命中计数列（列名 `{primer_id}_{LABEL}`）

## CLI 参考

```
Usage: seq_matcher [OPTIONS] --primer-csv <PRIMER_CSV> --library-csv <LIBRARY_CSV>

Options:
  -p, --primer-csv <PRIMER_CSV>              引物 CSV 文件路径
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
2. 构建 **Aho-Corasick 自动机**，将全部文库变体（原始 + 反向互补）编码为多模式匹配机，启动时构建一次，`Arc` 跨线程复用。
3. 目标序列按可配块大小分片，Rayon 多线程并行处理。
4. 每条序列**首次命中**引物即停止检索（first-match-wins）。
5. 匹配到引物后，对序列执行**单次 Aho-Corasick 扫描**即可检出所有命中变体（去重），替代逐变体 O(V) 次 `contains()` 调用。
6. 各线程内部无锁统计，块处理后合并至全局结果。
7. 实时进度条显示百分比、吞吐量、预计剩余时间。
8. 批量输出两个 CSV 文件。

### 性能

| 优化项 | 实现 |
|--------|------|
| 反向互补 | 编译期 `const` 128 字节 LUT，单周期映射 |
| 变体匹配 | Aho-Corasick 多模式搜索，O(L+M) 替代 O(V×L) |
| 并行处理 | Rayon work-stealing，chunk 级并行 |
| 线程间共享 | `Arc` 零拷贝传递引物、文库、自动机 |
| 输出 I/O | `BufWriter` 批量写入 |
| 进度反馈 | `\r` 原地刷新进度条，无额外 I/O 开销 |

## 依赖

| crate         | 用途                  |
|---------------|----------------------|
| clap          | CLI 参数解析           |
| csv           | CSV 读写              |
| rayon         | 数据并行               |
| serde         | 序列化 (derive)       |
| anyhow        | 错误处理               |
| aho-corasick  | 多模式子串匹配 (变体搜索) |

## License

MIT
