# smallRNA / PANDORA 分析模板

复制本目录的 `run.sh`、`shared/` 和 `config/` 到新课题，修改数据集配置和权威样本表后运行。现有课题使用自己的本地冻结代码；中心模板更新不会改变它们。

```text
smallrnaseq/
├── run.sh
├── shared/                 # 本项目冻结实现
└── config/                 # 单数据集；多数据集则放在 MII/config 等目录
```

## 使用

```bash
bash run.sh --list
bash run.sh --project . --stage describe --check
bash run.sh --project . --stage describe --modules qc,sncrna
bash run.sh --project MII --stage differential --modules mirna
```

入口按自身位置定位本地实现，命名数据集相对入口目录解析；`.` 表示入口目录自身。允许绝对数据集路径用于隔离验证。可从任意工作目录启动。`--list`、`--help` 和 `--check` 不写结果或缓存；预检通过仅表示输入可读取，不代表科学检验或结果已完成。

配置中的 `DB_ROOT`、`REFERENCE_ROOT`、`SPORTS_PL`、可选 `EXOGENOUS_DB_ROOT` 和 `ENDOSIRNA_REFS` 前缀均相对**具体数据集根**解析。样本表的 `raw_fastq_gz` 也遵循同一规则。不要填写用户 home 的固定路径；资源缺失应修配置，不建立旧目录跳转或回退。已有项目的集中参考位于 workspace 的 `reference/smallrna`，SPORTS 位于 `tools/sports1.1`。

## 配置与科学边界

`PROJECT_SPECIES` 必须显式为 `Mus musculus` 或 `Homo sapiens`，`RUN_ENDOSIRNA` 必须显式为 `0` 或 `1`；模式只由该字段确定。鼠使用 `SPORTS_LAYOUT=nested`，人使用 `flat`。参考映射格式为 `ENDOSIRNA_REFS='标签=相对FASTA前缀;标签2=相对前缀2'`，前缀不含 `.fa/.fasta`。

阶段为 `upstream、describe、differential、enrichment`。默认仅 describe；差异与富集需要明确选择且满足真实设计。语义模块为 `qc、pca、sncrna、tsrna、rsrna、mirna、pirna`；启用 endo 增加 `endo`；人类增加 `snrna`；人类描述步骤结果连续编号 `01_setup` 至 `09_snRNA_compare`，启用 endo 时使用独立 `10_endo_compare.Rmd` 与 `10_endosirna`，鼠原有步骤编号保持；人类方法比较从 `templates/integration/method_comparison` 复制到项目物种目录的联合分析中维护。单组学入口不再执行 `similarity、spikein、capture、capture_bias`。少于 3 个库不调度 PCA；人类目前没有 enrichment 阶段。

共同分类与原分析算法保持不变。endo 候选必须为 21–24 nt、宿主基因组匹配、参考全长零错配双向匹配，且原注释允许转类；显式 sncRNA 冲突保留原分类。候选资格不证明 Dicer/AGO2 依赖或功能。PN5 不默认剔除优势 18 nt 序列。SPORTS 分数分配和 reads 守恒继续审计；母类与子集不可重复累加。

serum 的供体 2 方法比较、供体 3 spike-in 工艺 QC 独立选择，均保留技术文库边界。目录归组不增加生物学重复数。单库仅支持描述；PCA 结构不等于独立验证。

## 输出与验证

新运行写入 `result/<阶段>`，包含 `run_status.tsv`、逐模块状态、输入与源码签名及 `index.html`。仅已选择且实际成功的模块写 completed。旧结果完整保存在各项目 `history/pre_layout_20261011`；旧签名不改写为新位置，历史 completed 不代表新入口已执行。

```bash
python3 scripts/make_fixture.py --output '/tmp/my smallrna fixture' --species Homo --endo 1
bash '/tmp/my smallrna fixture/run.sh' --project . --stage describe --modules sncrna,endo
```

合成生成器自建 SPORTS/参考/endo 表，不读取真实项目。它用于验证路径、互斥分类和冲突边界，不能证明原始 FASTQ 上游已经复现。执行结果与未跑边界见 [VALIDATION.md](docs/VALIDATION.md)。
