# SCI 路由与学习入口核验

核验日期：2026-10-14。状态：配置安装及本地验收通过。范围为 SCI 规则分工、学习导航、方法到绘图实现的连接；未执行真实项目分析或重新绘制科研图件。

## 当前实现

- 全局 `/home/h1028/.codex/agents/` 已更新19份文件并新增 `sci_core_contract.md`，安装后20份文件逐字节匹配审阅候选。14个角色的名称、描述、模型、推理强度及权限未变，reviewer仍为只读。
- 主SCI的 `developer_instructions` 从17,747字符减至4,246字符。共同科学原则集中在core；分析、图形、综述和排版规则按实际任务加载。主SCI管理路线和验收，详细规则留在对应角色与契约。
- 文章形式与研究类型分开。默认从当前问题选择一条主线与少量有明确用途的扩展，再核对输入资格；已有项目入口优先。
- 精选来源由evidence提炼，方法由匹配专科实现和验证，主SCI串行登记。学习库README是操作入口，`implementation_refs`连到代码，绘图目录的`method_ref`连回方法依据。
- 生物医学导航由正式卡元数据离线生成，共311项（263来源卡、23仓库、25综合）。32张来源卡缺少具体研究设计/原始类型，仍可检索并保留`unresolved`；全目录研究设计缺项共78项。没有刷新原阅读日期或将导航整理写成重新精读。
- 方法索引82个唯一锚点；绘图目录15条实现，12条已有示例、3条保持`source_audited / not_run`。预览导航从目录生成，核验既有图件回执，不执行统计或绘图重跑。

## 自动检查

以下命令从工作区根运行；配置检查先针对候选执行，再针对已安装全局配置执行。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s templates/tests -p test_sci_agent_config.py -q
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s sci_learning/tests -q
python3 sci_learning/scripts/build_catalogue.py --check
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s templates/figures/tests -q
python3 /home/h1028/.codex/skills/.system/skill-creator/scripts/quick_validate.py skills/omics-figures
```

| 检查 | 实际结果与范围 |
| --- | --- |
| SCI配置 | 19项通过；TOML、角色集合、权限、契约引用、领域路由与共同科学边界声明 |
| 学习目录生成器 | 11项通过；去重、分类轴独立、未知保留、只读检查、导航一致、启动目录无关及历史材料保护 |
| 学习导航一致性 | PASS，311项，`changed: []` |
| 图册生成器 | 7项通过；目录驱动、HTML转义、路径约束、回执与hash、未运行状态、重复生成一致及科研图件保护 |
| 绘图skill | `Skill is valid!` |
| 独立入口复核 | 82个方法锚点唯一，实现链接存在；绘图目录到方法依据的引用有效，职责无冲突 |

共37项自动测试通过。测试数量不包含skill格式检查、导航一致性检查及下列演练。静态配置检查只验证声明，不证明模型实际遵循规则。

## 六项独立上下文演练

任务及验收标准见[场景定义](../../tests/sci_routing_scenarios.md)。六次均使用fresh context，只提供假设任务和候选配置入口，不读取配置测试，不运行科研分析、联网或改文件。评估方式为当前会话默认模型读取候选指令，**并非配置中 `gpt-5.6-sol` 的生产调用或多轮性能测评**。本表的通过仅指对应场景的路线、实际读取范围和边界符合标准。

以下配置名均指当时临时候选，现已安装为全局同名文件；工作区路径保持当前入口。来源卡全文指卡片内容，不表示本轮重读原论文。

| 场景 | 实际读取的关键范围 | 观察结果 |
| --- | --- | --- |
| 母龄多时期RNA/smallRNA规划 | sci/core/article/analysis/analyst；学习README、相关索引；母龄综合、PMID 36697412与39367016卡；纯生信B04/B06 | PASS：一条阶段主线，年龄×时期及序列家族扩展；识别poly(A)测量层与非配对限制，拟议方法为`not_run`；未读图形生产/排版契约 |
| 每条件一个混池单细胞文库 | sci/core/article/analysis/singlecell；学习README及重复感知综合；绘图综合、目录与README；recipes.R的marker/program_matrix函数 | PASS：选择程序矩阵、捕获组成等描述图并给源表列；细胞不能补生物学n；标出WT/PFOS示例条件；未读排版契约 |
| GRCh38非配对多组学多轨图 | sci/core/analysis/figure/multiomics；方法综合、Signac卡；绘图目录、dense_recipes.py、输入检查与profile相关节 | PASS：实际发现GRCm39坐标与固定模态/配色限制；列出组装、区间、单位、方向与供者背景要求；不把坐标同位写成配对或因果；未读排版契约 |
| 公共转录组Research Letter | sci/core/article；学习README、纯生信与短文手册 | PASS：`computational_discovery`＋`research_letter`，从已有结果进入写作；列出待读的项目材料，不把结果存在当作验证；未加载完整分析/图形/排版契约 |
| 本地PMID 39367016学习 | sci/core/article/evidence；来源卡、既有学习记录、母龄综合相关节、方法索引、smallRNA模板README与B10路线 | PASS：提炼完整背景→序列家族→候选→独立验证；找到部分既有实现，区别已记录支持与未运行环节；明确本轮未读原论文；未加载分析/图形/排版契约 |
| 新特色方法沉淀 | sci/core/article/evidence；学习README、绘图README与纯生信B11相关节 | PASS：evidence提炼→专科最小可重放实现→主SCI串行登记；源码审阅/示例执行/项目验证分开，缺资源保留具体缺口与`not_run`；未加载无关分析/排版契约 |

原演练可按以下本地会话ID追溯，不复制整套会话日志进入工作区：

| 场景 | 会话ID |
| --- | --- |
| 母龄 | `01a135ff-c1cc-7140-8e81-56e8e9d4a864` |
| 单细胞 | `01a135ff-2270-75d3-86f4-fe5e1f4323a8` |
| 人源多轨 | `01a13600-b495-7680-a2a3-0c42f7687cbb` |
| Research Letter | `01a135fd-f60e-71e0-8bfd-88dd92279277` |
| 本地论文学习 | `01a13601-5931-7f73-b56b-0316f3582c4c` |
| 新方法沉淀 | `01a13602-2fad-7f80-8aa3-3973c3ed1265` |

## 核验边界与后续维护

这些检查支持配置可解析、入口一致，以及六个只读场景中按问题召回与限定推断。尚未量化真实SCI模型的长任务质量、耗时或上下文成本；主指令变短不等同于这些指标已经改善。示例绘图QA不等于上游论文复现或项目科学验证。

后续遇到实际错误时更新负责该规则的唯一入口，再复核对应场景；不重新堆叠主提示词或新增平行总索引。继续学习时更新原卡/手册和已有实现目录，由生成器维护派生导航。

## 安装内容指纹

以下SHA-256取自安装后的全局文件，供核对本次实际版本，不表示这些文件是运行结果。

| 文件 | SHA-256 |
| --- | --- |
| `sci.toml` | `f3cf7f90f47fc472f73ebb0f6557b036d4110b0b16b8b3c8c96efdf46598b9a1` |
| `sci_analysis_contract.md` | `830d7b2e97febd9c33f263a338a3fef452ec3171a07144443129eea65c836070` |
| `sci_analyst.toml` | `bdf7a81984400601adccbf82cede83ae380a2689de5288b7c3ce521ec2bf9c04` |
| `sci_core_contract.md` | `0654bbabda83d2f3aaa7f0b924231043f7c833cd41e0f14835865b2bde503c83` |
| `sci_epigenomics.toml` | `5f3549b727bee39db42c794b85eb0a7d385944f5d1dc2528b7f4f7f443ca16dc` |
| `sci_evidence.toml` | `e38d30863cf12f35dba55aafeb56a6a88fba9d7062f96fe546de3d2331ce96de` |
| `sci_figure_contract.md` | `980e66aadb4d8f3e30c92b786ef4c3c464f443da69cf3c4ff8da406d9284fff1` |
| `sci_manuscript_typography_contract.md` | `761836cd9a6ccfde539db9b45e616dbafad55a6f4d3112b0db19d438e534618d` |
| `sci_metabolomics.toml` | `d532c415049345baf2f02327271458db05ed77c8ffed99d87530ad7fe7cc544c` |
| `sci_ml.toml` | `8ff3ee49843c48750eae27971f3b775d9507acbe8f158d82061e62c85e29d7f0` |
| `sci_molecular_modeling.toml` | `9e971cf302e82aeca349be8d4b1c88d069d1fc2f8f8afb1e7a95fa19b5dad1d1` |
| `sci_multiomics.toml` | `1c04142bd033dfd3d1db3a43fee64c96f2dbe823538c7527ab9664092ccd1603` |
| `sci_network_pharmacology.toml` | `b161b59fb7fb36646b56375bd3aa83376d7f13392fdcd9940f1239bc9b3810b8` |
| `sci_review_figure_contract.md` | `918ba50aeb795f84cb63e5c03ed641a15a1b2bd004d2576e2121936e5aa3ed18` |
| `sci_review_writing_contract.md` | `ae1b0634b4028d4ec4eed7b5323650d40b0d5bb9879eeb80a6afbd44045d5510` |
| `sci_reviewer.toml` | `9a81714b6d51bc2a065a6115fdb060aef10ff9d06fef57769f228da1c379ade1` |
| `sci_singlecell.toml` | `6e9a3a290a35131405b0249fca7b3edfe354ef04948a06ff206503d50eff00a7` |
| `sci_spatialomics.toml` | `9e7c371b63dbc40570ef330b5eca13340666cf75f95e4e64f475c2a101ad2e11` |
| `sci_transcriptomics.toml` | `ad33ff922047d0669f95456f0dda234d5e7d1a51f9fc6a1f0d2cf85e7fba5bce` |
| `sci_writer.toml` | `372418f9707dc2ff46f9c4b4b4fabd178d9ccc567228c7bc4dca2187b78769cf` |
