# 可重放验证与留存证据

这些小脚本只生成合成数据或检查接口，源码路径由脚本自身位置确定。它们不依赖真实项目，也不读取既往临时目录；运行输出写入新建临时目录，不覆盖 `evidence/` 的本轮记录。

```bash
# 22 项输入、样本身份、行序、STALE、help/list/preflight、FASTQ 命名检查
python3 validation/contract_checks.py
# 基础 R 即可执行：rMATS 身份、两种 b1 方向、原始字段、NA 与反向 DE
Rscript validation/direction_checks.R
# 可选：使用现有 R 依赖重放 500×4 counts→DESeq2→01–05→HTML，约数分钟
bash validation/replay_core.sh
# Rscript 不在 PATH 时可显式给出可执行文件
RSCRIPT_BIN=/path/to/Rscript bash validation/replay_core.sh
```

`generate_counts_fixture.R` 保存了本轮随机种子 17、500 个 mouse Ensembl ID、4 个独立样本、60 个植入上调特征和交错样本表的生成规则。基因列表取自已安装 `org.Mm.eg.db`；包版本变化可能改变列表和结果，因此重放比较方向、样本范围、字段及状态，不要求报告字节完全一致。流程不会自动安装包。

`evidence/` 留存实际完成运行的 22 项边界结果、核心日志、模块执行状态、方向审计、分析参数、R 环境和小型汇总。日志中的临时路径是历史执行位置，不是运行依赖。`source_hashes_core_smoke.tsv` 绑定已完成核心运行的临时副本源码，不把后续轻量接口修订冒充为再次全量运行。没有复制大图或真实项目数据。`direction_checks.log` 是整理可重放脚本时对基础 R 小表检查的再次验证。

语法检查、合成小数据闭环和真实数据复现是不同证据等级。FASTQ 上游与 06–18 可选模块仍为 `not_run`；此前已验证的纯表报告情况见上级 `VALIDATION.md`。本目录不引入 R 测试框架。
