# serum 配对捕获与 spike-in 示例

这里保留“同供体成对分装属于工艺比较”“一个 spike-in 单点只能报告观察信号”的方法边界。它不进入通用独立两组流程，不把同供体分装算成独立生物学样本。

`paired_process_qc.R` 接受一个 TSV（pair_id、biological_sample_id、library_id、arm、clean_reads、spike_reads、spike_input_fmol），arm 必须为 control/test。输出每库 spike RPM、各对 test−control 的观察信号差以及单位清单，不产生人群 P 值；单序列单剂量信号不能换算分子回收率，也不能据此跨方法归一化。

```bash
Rscript paired_process_qc.R example_samples.tsv /tmp/serum_qc_example
```

example_samples.tsv 是人工合成输入，数值不来自真实研究。原始多模块实现的来源哈希记录在模板 provenance.tsv。完整 paired edgeR 技术捕获模型需要独立设计合同，本示例只提供可运行的描述性 QC。
