# PANDORA / SPORTS small RNA 模板

从 MII 项目提取现有 SPORTS 包装器和 9 个 R/Rmd 模块，保留互斥主类别、tRNA 嵌套类别与原始分母检查。复制本目录，修改 `config/project.env` 和 `config/samples.tsv` 后使用。样本不会按文件名推断组别、年龄或比较方向；没有 PN5 的特定序列删除逻辑。

## 配置与运行

`project.env` 是由入口 source 并导出的 Shell 参数文件，指定物种、组织/阶段、case/baseline、数据库前缀、接头、资源和各开关。所有相对路径相对于复制后的模板根目录。samples.tsv 每行一个文库，保留原 SPORTS 的 `sample/group` 命名；另给 `biological_sample_id/unit_type/batch` 定义独立单位。空 age_group 保持空，不猜 old/young。必须显式给出 group/condition，两个字段一致。

```bash
./sports.sh --help
./sports.sh --dry-run
./sports.sh --check
./sports.sh                    # 原始 FASTQ → cutadapt → SPORTS，可选 endo-siRNA
./run_rmd.sh --check            # 检查配置和显式 SPORTS 表输入，不写结果
./run_rmd.sh --run-id analysis1 # 默认模块 1 → 表达矩阵 → 描述图和中文报告
./run_rmd.sh --modules '2,3,4' --run-id analysis2
./run_rmd.sh --all --run-id all_modules
```

帮助、列表、dry-run 和 check 在创建结果目录前完成。REPORT_FIGURES=0 可生成仅表汇总报告。默认只跑最小闭环；2–8 为已有专题模块，9 由 RUN_ENDOSIRNA 决定。模块编号接受 08/09，不做八进制解析。每次下游运行使用新 run-id；已有 run-id 拒绝。原始 SPORTS 执行同样拒绝已有 SPORTS_RESULT_DIR，不因看到旧文件就声称当前配置完成，也不删除旧结果。设置新的输出路径重新执行。

原始 FASTQ 的 raw_fq 路径显式提供，文件名可与 sample 不同；不尝试旧文件名、后缀或目录回退。SPORTS_PL、参考 FASTA/Bowtie1 前缀、每类 small RNA 数据库、真实建库 ADAPTER3 必须配置。不同物种必须同时替换所有参考及 miRNA/OrgDb 参数；不能仅改 PROJECT_SPECIES。

## 统计边界

默认 RUN_INFERENCE=0，仅报告文库/组描述，不执行条件 DE。开启后仅接受每组至少 3 个已知独立 individual/pool、每单位一个文库、一个已知 batch 的独立两组设计。一个 pool 计一个单位；重复 technical library、未知身份和 batch 或批次模型需求均 NO-GO。本版不自动合并 SPORTS 注释后的技术库，也不支持配对/重复测量推断；需要对应设计时保持描述或另行制定方案。

实际 edgeR 使用显式 case − baseline，对表格重排不敏感；默认 |log2FC|≥log2(1.5) 且原始 P<0.05 为探索性筛选，保留 FDR。SPORTS 部分注释 Reads 是分配权重的小数，不能直接解释为独立原始计数；这些层会明确 skipped，**不四舍五入后检验**。整数序列/类别输入在其余资格通过时可运行 edgeR。资格和模块执行分别记录；GO 不等于已完成检验。

PCA 保留全部事前合格样本，用于 QC，不按组间分离程度选样本。主类别分母为 clean reads，并检查互斥类别回收 100%；tRNA 末端是嵌套子集，不重复相加。缺失原始 FASTQ 时原始读段 QC 为 missing_fastq，不能把 SPORTS 表输入闭环称为原始数据复现。

## 可选证据与示例

RUN_ENDOSIRNA 默认 0，无默认参考；打开必须给 ENDOSIRNA_REFS（上游 tag=Bowtie 前缀）及 ENDOSIRNA_REFERENCE_TSV（下游 ReferenceSet/ReferenceFile/UseInExactProxy/Reason）。exact sequence proxy 与 Bowtie reference mapping 分别保留证据级别，只支持候选，不证明内源来源或机制。所有参考来自显式配置，不搜索邻近项目。

RUN_TARGET_ENRICHMENT 默认 0，避免默认访问在线 miRNA 靶标服务。启用前设置物种、数据源和所需依赖；映射/预测/富集证据不能等同机制。MII 衰老 marker 文献比对和生物学解释移入 examples/mii，不混入通用 PCA 报告。

examples/serum 提供可独立运行的同供体成对工艺与 spike-in 观察信号 QC；不把分装当生物学重复、不把单剂量单序列信号换算为分子回收率。它与核心独立两组入口分开。

## 输出和验证

下游结果在 `ANALYSIS_ROOT/run-id`：`1_setup` 有注释和序列表达矩阵、来源和分母审计；`audit` 有冻结配置/样本、源码哈希、模块状态及推断资格；`sample_composition.csv`、`group_composition_descriptive.csv`、PDF/PNG 与 `report.html` 组成最小闭环。2–9 的细节输出沿用编号目录；未选择模块为 not_run。运行失败返回非零并标记 failed，不吞掉 Rmd 错误。

```bash
python3 scripts/make_fixture.py --output /tmp/smallrna_fixture --per-group 3
/tmp/smallrna_fixture/run_rmd.sh --run-id synthetic
Rscript examples/serum/paired_process_qc.R examples/serum/example_samples.tsv /tmp/serum_example
```

fixture 是独立复制目录，合成 SPORTS 三类表，不读取真实项目数据；默认拒绝覆盖。依赖说明见 DEPENDENCIES.md，实际结果与 not_run 边界见 VALIDATION.md，来源文件和哈希见 provenance.tsv。
