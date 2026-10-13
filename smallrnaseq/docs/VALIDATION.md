# smallRNA 当前职责与核验

新增物种层后，中心和三个课题的单组学实现继续独立冻结。当前29份入口、预检、上游和R/Rmd源码身份见 [SOURCE_SNAPSHOT.tsv](../SOURCE_SNAPSHOT.tsv)；方法比较四模块已抽离至中心integration/method_comparison及血清本地副本。

| 当前实际核验 | 结果 | 科学边界 |
|---|---|---|
| 六数据集双启动目录预检 | 12/12 PASS | 只读输入可用，不代表各阶段重跑 |
| 四物种/endo合成案例，在含空格新根执行 | 32/32 PASS | 分类、reads守恒、显式身份及签名；不是FASTQ复现 |
| MII miRNA差异 | completed；8/8图技术QA PASS | 31成熟miRNA表与前轮逐字节一致；MO−MY |
| MII富集入口 | completed；6/6图技术QA PASS | 当前DE签名available，但靶标0，GO/KEGG六项not_run |
| 血浆qc、sncrna | completed；42/42图技术QA PASS | 单库描述，联合资格not_applicable |
| 血清单组学全describe | 八模块completed；310/310图技术QA PASS | 12库描述，方法差异由独立联合入口负责 |

本轮真实记录集中于workspace的 [maintenance/species_organization/smallrna_sc](../../../maintenance/species_organization/smallrna_sc/README.md)。MII实际producer门验证了current通过、failed拒绝、MD5变化拒绝和缺输入拒绝。没有安装或升级依赖，没有重跑原始FASTQ、SPORTS或长期上游。

本目录validation下的旧运行记录保持原始身份；[前轮验收文档](../history/pre_species_organization_20261011/VALIDATION.md)描述抽离前33份源码及原血清差异入口，不能作为现役单组学仍执行方法比较的证明。其他三个母龄时期只做预检，历史describe凭据不冒充当前源码completed。实际运行与技术QA不等于生物学结论或人工审美验收。

2026-10-12 编号调整核验：人类主步骤结果为 01–09，可选 endo 为 10；鼠原编号保留。四个 Mus/Homo × endo0/1 隔离案例重新实际执行，32 项检查通过，116 个复制源码与当前模板逐字节一致；包括含空格路径、只读预检、缺失 FASTQ 拒绝、无效模块拒绝、reads 守恒及 sncRNA 冲突边界。实际记录见 [template_portability.json](../../../maintenance/analysis_completion_20261012/smallrna/template_portability.json)。
