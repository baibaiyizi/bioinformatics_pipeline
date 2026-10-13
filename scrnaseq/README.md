# scRNA-seq 独立两组模板

复制本目录后，配置 `config/project.yml` 与 `config/samples.tsv` 即可从原始文库矩阵或 FASTQ 开始。核心支持 Cell Ranger → QC → Scanpy PCA/UMAP/Leiden → 可配置 marker 注释 → 文库描述 → 可选样本级 pseudobulk/PyDESeq2 → 中文 HTML 报告与执行审计。默认不校正批次、不挑选能分开条件的样本、不启用条件推断。

绘图使用本地 `src/scrna_template/figure_style.json` 冻结配置，通用入口见[生信绘图库](../figures/README.md)。对照/处理名称从现有statistics.control/case绑定蓝/粉，文库色从样本表condition绑定，不依赖文库名。cell type和cluster按profile固定名称；新类别须显式登记，不能按出现顺序重配。PFOS专题可选源码在 `examples/islet/optional/`，其真实107图验证不能冒充通用模板全流程复现。

## 环境和入口

核心 Python 要求 `>=3.13,<3.14`，依赖权威为 `pyproject.toml`；`requirements-validated.txt` 记录本次环境的精确版本。可在自己的环境中安装 `pip install -e '.[replicated,test]'`，仅描述时省略 replicated。服务器原环境未升级或安装软件。

```bash
# SCRNA_PYTHON 可设为已准备好的 Python 3.13 绝对路径；不设置则用 PATH 中 python。
./run_scrna.sh --help
./run_scrna.sh --list-steps
./run_scrna.sh --config config/project.yml validate
./run_scrna.sh --config config/project.yml validate --full-input-audit
./run_scrna.sh --config config/project.yml run --run-id analysis_001
./run_scrna.sh --config config/project.yml plots --source-run runs/analysis_001 --check
./run_scrna.sh --config config/project.yml plots --source-run runs/analysis_001
```

配置路径相对于显式 `--config` 调用路径解析；`project.root: ..` 相对于 YAML 所在目录解析，其余数据、参考和结果路径相对于 project.root。帮助/步骤列表不读取配置；validate 为只读，不创建结果。每次分析必须使用新 run-id，已存在目录会拒绝，避免把旧文件当成当前输出。

## 样本表和统计单位

每行一个技术文库，必备列为 `library_id, biological_sample_id, unit_type, condition, sample_batch, chemistry, matrix, r1, r2, i1`。TSV 的可选空字段保持为空；身份未知填 unknown，绝不推测分组或动物身份。

- `library_id` 必须唯一且为安全文件名。`biological_sample_id` 是独立动物/人或一个独立混池的 ID。同一生物学样本多个技术文库使用相同 ID；`unit_type` 为 individual、pool 或 unknown。每个 pool 只计一个样本，不使用池内动物数扩增 n。填写单位即声明样本独立；家系/窝效应等聚类设计不在本版支持范围。
- `condition` 必须恰好对应配置中的 case 与 control，方向恒为 case − control，不受表格顺序影响。`sample_batch` 是样本级、可评估的批次因素，同一样本多库必须一致；若技术文库跨批次，本版不支持其复杂结构，降为描述。
- 只有独立两组、每组至少 3 个已知单位、单位/条件/批次一致、设计满秩并有残差自由度时允许推断。资格同时在全设计和每细胞类型最低细胞数过滤之后检查。未知单位/缺失 batch/样本不足/设计混杂仅做描述；显式 paired/repeated 设计及非空 subject_id、pair_id、timepoint 字段直接拒绝。
- 先按 biological_sample_id × cell_type 汇总全部技术文库的 `counts_raw` 非负整数，再使用 `~condition` 或 `~sample_batch + condition`。不对归一化表达做 DESeq2，不把细胞作为 n。

## 两类显式输入

`inputs.mode: h5ad`：每个 matrix 指向独立文库的原始 counts AnnData。若有 counts_raw 层则读取该层，否则读取 X 并验证非负整数。行是细胞、列是基因 symbol，索引唯一；文库须使用相同参考和完整基因集合。不会读取来源项目现成 atlas/归一化结果作为原始计数。matrix 模式记录 SHA-256、完整计数和索引检查，不声称通过 FASTQ MD5 审计。

`inputs.mode: fastq`：r1/r2 列显式列出原始 gzip FASTQ；同文库多 chunk 用分号分隔，R1/R2 数量必须相同。i1 可空，有则数量与主读段相同。供应商文件名可以不同于 library_id。MD5 清单必须恰好覆盖所列 FASTQ；basename 重复拒绝。正式运行检查最后修改年龄、无活动写句柄、完整 MD5、全长记录/序列质量长度/R1-R2（及可选 I1）header 配对，并检查校验前后大小和 mtime 不变。FASTQ 原件不改名；审计后只在本 run 下创建规范名字的符号链接供 Cell Ranger 读取。无写句柄检查需要 lsof，不能完成时明确失败。

Cell Ranger 可执行程序和转录组参考由配置给出；资源和每文库 chemistry 显式配置。FASTQ 上游在本轮未跑真实数据，状态 `not_run`，不会因矩阵入口通过就宣称完成 Cell Ranger 复现。

## 分析和结果

QC 根据事前配置的 UMI、检测基因数、线粒体比例过滤。Scrublet 可显式启用并给出建库对应 expected_doublet_rate，失败终止；默认关闭。环境 RNA 校正本版未执行，报告明确记录原始计数保留。marker 文件提供 positive/negative 集，使用来源实现的 cluster mean log1p 基因标准化证据；低分差/低覆盖为 unknown，空 marker 配置全部 unknown。未知标签不进入细胞类型 DE。

`statistics.inference_enabled: true` 才运行 PyDESeq2。默认 `P < 0.05` 且 `|log2FC| >= 0.58` 为探索性筛选，完整表同时保留 pvalue/padj 和 raw-P/BH 两个标志；NA 不改为 0。cluster marker P 值仅作注释证据，不能解释为条件差异。组成输出为捕获组成描述，不进行总体比例推断。

主要输出：`00_audit` 保存输入哈希、资格、阶段状态、运行结果和输出哈希；`05_descriptive` 保存文库组成、样本级整数聚合和描述 CPM 比值；`06_optional/inference_status.tsv` 分开记录 eligibility 与 execution（disabled/skipped_ineligible/completed/failed）；实际完成才存在 DE 表。某细胞类型推断失败会保留失败状态，生成可查看的报告，并使整个命令非零退出。`09_reports/scrna_report_zh.html` 支持有图、仅表、跳过模块与失败记录。

纯绘图入口从原运行 `output_sha256.tsv` 核验每张图的数值源，并要求 producer 为 PASS；只支持通用核心真实声明的 UMAP 与捕获组成图，遇到专题扩展图型直接拒绝。它不重跑统计，也不写入原 `runs/`。发布的 `result/` 按 `00_audit`、`01_qc`、`02_atlas`、`03_annotation`、`04_descriptive` 分层，每层分别保存 `plots/`、`tables/`，总入口为 `result/index.html`。未产生的阶段不建空图表。

核心绘图复用本地 `figure_style.py`：183 mm、Liberation Sans、8/7 pt、300 dpi RGB 白底 PNG 与嵌入字体 PDF。WT/PFOS 有明确语义配色，其余配置组沿用类别色；细胞类型配色由实际项目设置。专题示例的扩展纯绘图源码仍是参考档案，不能据此声称通用核心支持全部专题模块。

## 合成验证

```bash
./run_scrna.sh fixture --output /tmp/my_scrna_fixture --per-group 3 --technical-libraries 2
./run_scrna.sh --config /tmp/my_scrna_fixture/config/project.yml run --run-id synthetic
python -m pytest tests
```

生成器拒绝覆盖已有目录，内部构建配置，不读真实项目样本或配置。合成 Gene000 在 treated 中植入上调，只用于验证方向/计数/执行合同。`examples/islet` 保存胰岛专题参数和扩展源码；它们不作为已验证通用模块开放。实际测试结果和上游限制见[内部验证记录](docs/VALIDATION.md)，来源及 SHA-256 见 provenance.tsv。
