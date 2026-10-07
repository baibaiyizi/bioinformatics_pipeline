# 合成 consensus counts 示例

500 个模拟 BED 区间，4 个独立样本，整数 counts；前 30 峰 treated 放大 8 倍，31–60 峰 control 放大 8 倍。固定已生成矩阵用于检验方向与运行入口，不是生物学数据；真实样本效应/功效不可由此推断。样本使用 01–04，专门覆盖身份前导零。

复制整个模板至一个新的临时目录后，在 config/project.env 中设置以下三项，再运行：

```bash
ATAC_INPUT_MODE=counts
ATAC_SAMPLES=examples/consensus_counts/samples.tsv
ATAC_COUNTS_FILE=examples/consensus_counts/counts.tsv
# direction defaults: treated - control
bash run_rmd.sh --check
bash run_rmd.sh
```

上述前三项应写入参数文件，或以同名环境变量传入每次入口调用。预期 500 个检验峰；原始 p<0.05 主集合 72 峰（DESeq2 1.50.2）。前 30 个 Fold>1，31–60 个 Fold<−1。修改比较方向应在另一个新目录运行，避免把旧方向缓存纳入；所有 Fold 反向。过滤、归一化与 BH 结果来自实际 DESeq2。
